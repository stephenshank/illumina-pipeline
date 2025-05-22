sra_examples=("cambodia" "stephen-sarscov2")

for example in "${sra_examples[@]}"; do
  mkdir -p $example/data
  for accession in $(cat $example/ids.txt); do
    fasterq-dump -O $example/data $accession
  done
done