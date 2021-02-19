#! /bin/sh

#usare: ./costruisci_input.sh "IRC_OUTPUT.log" "number_of_atoms_involved"
#starting point: first point analyzed on PES 
#steps: number of point considered on PES

Atoms=$2
let "FCread=0"
let "FCreadmax = $Atoms*3/5"
for i in `seq 0 1 $FCreadmax`
do
let "FCread=$FCread+$Atoms*3-5*$i+1"
done
#echo $FCread

let "FCreadm1=$FCread-1"

let "Gread=$Atoms+2"
let "Sread=$Atoms+4"
echo $FCread
echo $Gread
echo $Sread

rm -f grad.txt
rm -f grad1.txt
rm -f grad2.txt
rm -f grada.txt
egrep -A$Gread "Forces " $1 > grad1.txt
tail -$Atoms grad1.txt > grad2.txt
echo "gradient" >> grada.txt
cat grada.txt grad2.txt >> grad.txt
rm -f grad1.txt
rm -f grad2.txt
rm -f grada.txt


rm -f force_constants2.txt 
rm -f force_constants1.txt 
rm -f force_constants.txt 
rm -f hessa.txt
egrep -A$FCread "Force constants in Cartesian coordinates" $1 > force_constants1.txt 
tail -$FCread force_constants1.txt  > force_constants2.txt 
sed -ie 's/D/E/g' force_constants2.txt
echo "Hessian " >> hessa.txt  
cat hessa.txt force_constants2.txt >> force_constants.txt
rm -f force_constants1.txt 
rm -f force_constants2.txt 
rm -f force_constants2.txte 
rm -f hessa.txt



rm -f geometries1.txt
rm -f geom.txt
rm -f geoma.txt
rm -f geomb.txt
egrep -A$Sread "Input orientation" $1 > geometries1.txt
echo "Step 1 " >> geoma.txt  
echo "geometry " >> geoma.txt  
tail -$Atoms geometries1.txt > geomb.txt
cat geoma.txt geomb.txt >> geom.txt
rm -f geoma.txt
rm -f geomb.txt
rm -f geometries1.txt

#sed -ie  '/ATOMS/i\'$Atoms'\' hind_rot_head.dat
#sed -ie  's/ATOMS/'$Atoms'/g' hind_rot_head.dat

echo 'Number_of_Atoms:        ' $Atoms > temp.log

cat temp.log hind_rot_head.dat hrdata4proj.dat geom.txt grad.txt force_constants.txt > RPHt_input_data.dat

rm -f geom.txt
rm -f grad.txt
rm -f force_constants.txt
rm -f temp.log

#egrep "SCF Done" $1 > Energies_sym.txt



