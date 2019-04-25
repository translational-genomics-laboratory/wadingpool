

for id in `cat $IDLIST | sort | uniq`; do
  DUPM=`ls $BAMDIR/$id*.dup_metric.txt`
  if [[ -f $DUPM ]]; then
    HEADER=`cat $DUPM | tail -n +7 | head -1`
    METRIC=`cat $DUPM | tail -n +8 | head -1`
    echo -e "ID\\t$HEADER" > $QC/$id.dup_metric.txt
    echo -e "$id\\t$METRIC" >> $QC/$id.dup_metric.txt
  fi
  PCI=`ls $PCQCDIR/$id.isize.txt`
  if [[ -f $PCI ]]; then
    HEADER=`cat $PCI | tail -n +7 | head -1`
    METRIC=`cat $PCI | tail -n +8 | head -1`
    echo -e "ID\\t$HEADER" > $QC/$id.insertsize.txt
    echo -e "$id\\t$METRIC" >> $QC/$id.insertsize.txt
  fi
  sWGSM=`ls $PCQCDIR/$id.wgsMetrics.txt`
  if [[ -f $sWGSM ]]; then
    HEADER=`cat $sWGSM | tail -n +7 | head -1`
    METRIC=`cat $sWGSM | tail -n +8 | head -1`
    echo -e "ID\\t$HEADER" > $QC/$id.sWGSMetrics.txt
    echo -e "$id\\t$METRIC" >> $QC/$id.sWGSMetrics.txt
  fi
done
TEL=$TL/allTelbamLengths.csv
cat $TEL | tr "," "\t" > $QC/telomereLengthMetrics.txt
