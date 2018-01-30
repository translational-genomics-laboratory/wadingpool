PDIR=$1
cd ${PDIR}/output_wgs
rm allOL.wgsMetrics.txt

for i in $(ls . | grep "wgsMetrics.txt"); do 
  head -7 $i | tail -1 > tmp.txt
  sed "s/^/SAMPLES\\t/" tmp.txt > header.txt

  head -8 $i | tail -1 >> wgsMetrics.txt; 
  echo $i >> samples.txt
done

paste -d "\t" samples.txt  wgsMetrics.txt > tmp.txt
cat header.txt tmp.txt > ../output_summary/allOL.wgsMetrics.txt

rm wgsMetrics.txt
rm samples.txt
rm header.txt
rm tmp.txt