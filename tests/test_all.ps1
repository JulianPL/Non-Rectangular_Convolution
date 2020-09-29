<# 
    Powershell-script for checking format and for unit-testing.
#> 


$Diff = yapf -d $PSScriptRoot"\..\nrconv\primes.py"
If ($Diff -ne $null)
{
    echo "WARNING: nrconv\primes.py not formatted correctly"
}

$Diff = yapf -d $PSScriptRoot"\..\nrconv\convolution.py"
If ($Diff -ne $null)
{
    echo "WARNING: nrconv\convolution.py not formatted correctly"
}

$Diff = yapf -d $PSScriptRoot"\..\tests\test_primes.py"
If ($Diff -ne $null)
{
    echo "WARNING: tests\test_primes.py not formatted correctly"
}

$Diff = yapf -d $PSScriptRoot"\..\tests\test_convolution.py"
If ($Diff -ne $null)
{
    echo "WARNING: tests\test_convolution.py not formatted correctly"
}

python3 -m unittest tests.test_primes 2> $null
If ($LASTEXITCODE -ne "0")
{
    echo "ERROR: tests\test_primes.py found an error:"
    python3 -m unittest tests.test_primes
}

python3 -m unittest tests.test_convolution 2> $null
If ($LASTEXITCODE -ne "0")
{
    echo "ERROR: tests\test_convolution.py found an error:"
    python3 -m unittest tests.test_convolution
}
