#!/bin/bash

for f in $( ls *.fasta);
do sed -i "s/\\./_/g" $f
   sed -i "s/~/-/g" $f
   sed -i "s/|/_/g" $f;
done

# ramxl 8.2.11 

for f in $(ls *.fasta);

do ~/yaofile/raxmlfile/raxml -f a -p 123 -s $f -x 666 -#autoMRE -m GTRGAMMA -n $f;

        for k in $(ls RAxML_bestTree.$f);
        do cat $f $k > ${f}t;
        done

mv RAxML_bipartitions.$f raxml_${f/.fasta/.tre}


done

rm RAxML*
rm *.fasta

filedir=/home/yao/yaofile/dS/20190123-n1

for f in `ls $filedir/*.fastat`;
do echo `( echo "1"; echo "3"; echo $f; echo "Y") | HYPHYMP /home/yao/yaofile/dS/script/dSdN.bf.txt` > ${f/.fastat/.dNdS} 
done 

rm messages.log

