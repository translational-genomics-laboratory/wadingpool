PDIR=$1
cd ${PDIR}/output_iSize
rm allOL.isize.txt

for i in $(ls . | grep "isize.txt"); do 
  head -7 $i | tail -1 > tmp.txt
  sed "s/^/SAMPLES\\t/" tmp.txt > header.txt

  head -8 $i | tail -1 >> isize.txt; 
  echo $i >> samples.txt
done

paste -d "\t" samples.txt  isize.txt > tmp.txt
cat header.txt tmp.txt > ../output_summary/allOL.isize.txt

rm isize.txt
rm samples.txt
rm header.txt
rm tmp.txt