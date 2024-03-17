#! /bin/bash
echo "run the COLONY analysis for 2010 to 2014"

for i in {2010..2014}
do
echo "running $i ..."
	./colony2s.ifort.out OFN:COLONY_output_"$i".txt IFN:COLONY_input_"$i".txt
done

