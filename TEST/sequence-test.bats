#!./TEST/libs/bats/bin/bats

load 'libs/bats-support/load'
load 'libs/bats-assert/load'

@test "Runs Test.fastq.gz" {
    ln -s TEST/Test*.fastq.gz .
    run ./bactofidia.sh Test
    [ "$status" -eq 0 ]
    rm Test*fastq.gz
}
