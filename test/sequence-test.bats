#!./TEST/libs/bats/bin/bats

load 'libs/bats-support/load'
load 'libs/bats-assert/load'

@test "Checks if bactofidia runs with test sample" {
    cp test/Test*.fastq.gz .
    run ./bactofidia.sh Test
    [ "$status" -eq 0 ]
    rm Test*fastq.gz
}
