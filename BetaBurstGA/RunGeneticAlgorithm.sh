#!/bin/bash
BASEDIR=$(pwd)
JPATH="$BASEDIR/Julia_Code/Scripts"
MPATH="$BASEDIR/HMM"
JULIA=julia

echo "Running inital Generation"

# Initial Population

$JULIA $JPATH/makeInitialGeneration.jl
$JULIA $JPATH/RunMassModel.jl
matlab -nodisplay -nojvm -nosplash -nodesktop -r "try;run('$MPATH/Run.m');catch;end;exit"

j=0
for i in {1..100}
do
echo "Generation $j: " >> bestFitnessLog.txt
echo Running $i th Generation
let j=$j+1
$JULIA $JPATH/makeNextGeneration.jl
$JULIA $JPATH/RunMassModel.jl
cat $JPATH/../tmpfiles/tmpBestFitnessLog.txt >> bestFitnessLog.txt
matlab -nojvm -nodesktop -nodisplay  -r "try;run('$MPATH/Run.m');catch;end;exit"
done
echo "Generation $numGenerations: " >> bestFitnessLog.txt
$JULIA $JPATH/comparefitness.jl
cat $JPATH/../tmpfiles/bestFitnessLog.txt >> bestFitnessLog.txt
rm $JPATH/../tmpfiles/bestTotFitnessPREV.jld

echo done.
