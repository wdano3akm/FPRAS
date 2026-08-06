#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: compile_complete_nnf.sh INPUT.cnf OUTPUT.nnf [D4_OPTIONS...]

Runs d4's d-DNNF compiler, then adds a tautology (x OR -x) for every
variable declared by the CNF but absent from the dumped NNF.

Set D4_COMPILER to select the compiler executable. It defaults to the
d4compiler executable in the repository root.
EOF
}

if (( $# < 2 )); then
    usage
    exit 2
fi

input=$1
output=$2
shift 2

if [[ ! -f $input ]]; then
    echo "error: CNF input does not exist: $input" >&2
    exit 1
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
compiler=${D4_COMPILER:-"$repo_root/d4compiler"}

if [[ ! -x $compiler ]]; then
    echo "error: d4 compiler is not executable: $compiler" >&2
    echo "set D4_COMPILER=/path/to/compiler to select it" >&2
    exit 1
fi

cnf_metadata=$(
    awk '
        $1 == "p" && $2 == "cnf" {
            if (found) {
                print "error: CNF contains more than one p cnf header" > "/dev/stderr"
                exit 1
            }
            if ($3 !~ /^[0-9]+$/) {
                print "error: invalid variable count in p cnf header" > "/dev/stderr"
                exit 1
            }
            count = $3
            found = 1
            next
        }

        $1 !~ /^c$/ && found {
            for (i = 1; i <= NF; ++i) {
                literal = $i + 0
                if (literal == 0) {
                    continue
                }
                atom = literal < 0 ? -literal : literal
                if (atom < 1 || atom > count) {
                    printf "error: CNF literal %d is outside its declared variable universe 1..%d\n", \
                        literal, count > "/dev/stderr"
                    invalid = 1
                }
                used[atom] = 1
            }
        }

        END {
            if (invalid) {
                exit 1
            } else if (!found) {
                print "error: CNF has no p cnf header" > "/dev/stderr"
                exit 1
            }

            printf "%d", count
            for (atom = 1; atom <= count; ++atom) {
                if (!(atom in used)) {
                    printf " %d", atom
                }
            }
            printf "\n"
        }
    ' "$input"
)

read -r -a cnf_fields <<< "$cnf_metadata"
declared_variables=${cnf_fields[0]}
cnf_free_variables=("${cnf_fields[@]:1}")
cnf_free_list=${cnf_free_variables[*]}

output_dir=$(dirname -- "$output")
if [[ ! -d $output_dir ]]; then
    echo "error: output directory does not exist: $output_dir" >&2
    exit 1
fi

raw_nnf=$(mktemp "$output_dir/.d4-raw.XXXXXX")
complete_nnf=$(mktemp "$output_dir/.d4-complete.XXXXXX")
cleanup() {
    rm -f -- "$raw_nnf" "$complete_nnf"
}
trap cleanup EXIT

"$compiler" \
    -i "$input" \
    --dump-file "$raw_nnf" \
    "$@"

if [[ ! -s $raw_nnf ]]; then
    echo "error: d4 did not produce an NNF dump" >&2
    exit 1
fi

# d4 numbers and prints the root before recursively printing its descendants.
# Record that first node, the largest node id, and all atoms present on arcs.
metadata=$(
    awk -v declared="$declared_variables" -v free_list="$cnf_free_list" '
        BEGIN {
            count = split(free_list, free_atoms, " ")
            for (i = 1; i <= count; ++i) {
                if (free_atoms[i] != "") {
                    free_variable[free_atoms[i]] = 1
                }
            }
        }

        $1 ~ /^[oatf]$/ && NF == 3 && $3 == 0 {
            id = $2 + 0
            if (!have_root) {
                root = id
                have_root = 1
            }
            if (id > max_id) {
                max_id = id
            }
            next
        }

        $1 ~ /^[0-9]+$/ {
            for (i = 3; i <= NF && $i != 0; ++i) {
                literal = $i + 0
                atom = literal < 0 ? -literal : literal
                if (atom < 1 || atom > declared) {
                    printf "error: NNF literal %d is outside the CNF variable universe 1..%d\n", \
                        literal, declared > "/dev/stderr"
                    invalid = 1
                }
                present[atom] = 1
            }
        }

        END {
            if (invalid) {
                exit 1
            }
            if (!have_root) {
                print "error: NNF contains no node declaration" > "/dev/stderr"
                exit 1
            }

            for (atom = 1; atom <= declared; ++atom) {
                if (!(atom in present) && !(atom in free_variable)) {
                    printf "error: d4 omitted constrained CNF variable %d; refusing to treat it as free\n", \
                        atom > "/dev/stderr"
                    unsafe = 1
                }
            }
            if (unsafe) {
                exit 1
            }

            printf "%d %d", root, max_id
            for (atom = 1; atom <= declared; ++atom) {
                if (!(atom in present) && atom in free_variable) {
                    printf " %d", atom
                }
            }
            printf "\n"
        }
    ' "$raw_nnf"
)

read -r -a fields <<< "$metadata"
original_root=${fields[0]}
max_node_id=${fields[1]}
missing_variables=("${fields[@]:2}")

if (( ${#missing_variables[@]} == 0 )); then
    mv -- "$raw_nnf" "$output"
    echo "wrote $output; all $declared_variables CNF variables were already represented" >&2
    exit 0
fi

required_nodes=$((2 + ${#missing_variables[@]}))
if (( max_node_id > 4294967294 - required_nodes )); then
    echo "error: adding tautologies would overflow the NNF node-id range" >&2
    exit 1
fi

new_root=$((max_node_id + 1))
true_node=$((max_node_id + 2))
next_node=$((max_node_id + 3))
tautology_nodes=()

for _ in "${missing_variables[@]}"; do
    tautology_nodes+=("$next_node")
    next_node=$((next_node + 1))
done

{
    echo "a $new_root 0"
    echo "t $true_node 0"
    for node in "${tautology_nodes[@]}"; do
        echo "o $node 0"
    done

    cat -- "$raw_nnf"

    echo "$new_root $original_root 0"
    for index in "${!missing_variables[@]}"; do
        variable=${missing_variables[$index]}
        node=${tautology_nodes[$index]}
        echo "$new_root $node 0"
        echo "$node $true_node $variable 0"
        echo "$node $true_node -$variable 0"
    done
} > "$complete_nnf"

mv -- "$complete_nnf" "$output"
echo "wrote $output; added ${#missing_variables[@]} free-variable tautologies:" >&2
printf '  %s\n' "${missing_variables[*]}" >&2
