# blast鉴定耐药基因三步曲
## 数据准备

-. 新下载的NCBI细菌数据库更新到服务器上
*/share/disk5/zhuqh/bacteria_ncbi_1_10_2019/database*，
大小为563G，包含95种菌，各level的strain数目统计结果在
*/share/disk5/zhengx/identify_ARGs_1-14/0base_file/122species_count.txt*
在这里先做那8种菌的耐药基因鉴定，八种菌名在
*/share/disk5/zhengx/identify_ARGs_1-14/0base_file/8species.name*
里。
- 首先将NCBI database里面的protein 序列提取出来，然后解压，放在
*/share/disk5/zhengx/identify_ARGs_1-14/1fasta_file.pro*，
在1fasta_file.pro下每种菌一个文件夹，数据结构为:**1fasta_file.pro/species.faa/level_GCA_species.faa。**
提取序列：               
`nohup cat /share/disk5/zhengx/identify_ARGs_1-14/0base_file/8species.name  | while read i ; do ls /share/disk5/zhuqh/bacteria_ncbi_1_10_2019/database/${i}/Complete/| while read j ; do cp /share/disk5/zhuqh/bacteria_ncbi_1_10_2019/database/${i}/Complete/${j}/*_protein.faa.gz ./${i}.faa/comp_${j}_protein.faa.gz  >>./nopro_fa.txt  2>&1`

(在这里，追加重定向到文件./nopro_fa.txt ，是为了将错误输出定向到文件中，因为database里有的stain没有相应的protein文件，因此将错误输出到文件中，方便统计有哪些strain没有protein文件） ; done;done &

- 统计提取出的序列和提取不出来的序列加起来数目对不对：

   - 统计已经提取出来的level 为chromosome的strain有多少
`cat ../0base_file/8species.name | while read i ; do a=`ls ./$i.faa | grep -c  "chro"` ; echo $i   $a >>1chro_count;done`
   - 统计报错信息中的没有成功提取的序列的数目，加起来看和总数对不对得上
`cut -d "/" -f 7 nopro_chro_file | sort | uniq -c | awk '{print $2 "\t" $1}' >2no-chro_count`

- **当提取第二部分序列（除了重要的八种菌之外的其他87种菌，有一个重大问题是统计文件数目，有没有全部cp完成）**
（wd：/p200/husn_group/zhengxin/identify_ARGs_1-14/count）

`nohup cat 3-7_bsi_all_species.txt | while read i ;do ls /p200/husn_group/zhuqh/bacteria_ncbi_1_10_2019/database/${i} |while read j ;do ls /p200/husn_group/zhuqh/bacteria_ncbi_1_10_2019/database/${i}/${j} | while read z ;do a=`ls /p200/husn_group/zhuqh/bacteria_ncbi_1_10_2019/database/${i}/${j}/${z} | grep -c "protein.faa.gz"` ;if [ $a -gt 0 ] ;then echo $z >>2all_strain_has_pro_id ;else echo $z >>2all_strain_nohas_pro_id ;fi ; 2>err.log; done ;done  ;done &`

      - 解压:`gunzip *`

## blastp鉴定耐药基因

`nohup ls ./import-temp_class1_species_1fasta_file.pro | while read i ; do name=$(basename $i .faa) ; ls import-temp_class1_species_1fasta_file.pro/$i | while read j ; do strain=$(basename $j _protein.faa) ; /p200/husn_group/zhengxin/software/ncbi-blast-2.2.28+/bin/blastp -query /p200/husn_group/zhengxin/identify_ARGs_1-14/import-temp_class1_species_1fasta_file.pro/${i}/${j} -db /p200/husn_group/zhengxin/CARD_seq/new-symbol_CARD_protein_new.fasta -out 2blast_out.e-20/${name}.blout/${strain}_pro.blout  -evalue 1e-20 -outfmt 0  -num_threads 8  ;done;done >>class1.blast.log 2>&1  &  `

#【鉴定耐药基因的参数，分为两步来卡，一步是在blastp，注意blastp的版本要是2.2.28的，因为Eblastx版本为4.1， 是根据2.2.28来做的，在blastp这一步只卡evalue，为1e-20，这里的参数是按照NCBI的耐药细菌数据库来卡的】

#【另外的参数identity为90%，qcoverage为80%，tcoverage为80%，这些都是在EblastX这一步卡的】

## Eblastx   

#【以Acinetobacter_baumannii为例，参数identity为90%，qcoverage为80%，tcoverage为80%】
`ls 2blast_out.e-20/Acinetobacter_baumannii.blout/ | while read i ; do perl /share/disk5/zhengx/software/EblastX_4.1.pl -i 2blast_out.e-20/Acinetobacter_baumannii.blout/$i -o 3Ebl_out/Acinetobacter_baumannii_Eblax/${i}.ebl -id 90 -qp 0.8 -tp 0.8 -b 1 -p 20  >>ebl_AB_out.log 2>&1 ;done`
