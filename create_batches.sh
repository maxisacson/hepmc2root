#!/bin/bash

create_pbs_job() {
    local g4mac=$1
    local output=$2
    local jobfile=${g4mac##g4_}
    jobfile=${jobfile%.mac}
    jobfile="pbs_${jobfile}.sh"

cat > "$jobfile" <<EOF
#!/bin/bash

#PBS -N ${g4mac%.mac}
#PBS -l nodes=1:ppn=1

cd \$PBS_O_WORKDIR
./GenericDetector $g4mac $output
EOF
}

input="$1"
base_output=${input##*/}
events_per_file=50

build/split_hepmc2 "$input" "$base_output" -e "$events_per_file" -n 200

while IFS= read -r -d '' x; do
    base=${x##*/}
    g4mac="g4_${base}.mac"

    echo "Creating job for $x"

cat > "$g4mac" <<EOF
/control/verbose 2
/run/verbose 2
/generation/select hepmcAscii
/generator/hepmcAscii/open ${x}
/generator/hepmcAscii/verbose 0
/run/beamOn ${events_per_file}
EOF

    create_pbs_job "$g4mac" "${base}.root"
done < <(find $(pwd) -type f -regex "^.*${base_output}\.[0-9]+$" -print0 | sort -z)
