#!./TEST/libs/bats/bin/bats

load 'libs/bats-support/load'
load 'libs/bats-assert/load'

@test "Runs Test.fastq.gz" {
    run ./bactofidia.sh TEST/Test
}
