import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings('ignore')
i = 0
with open ('/home/olesjabogomaz/grid_nodup_mm10_ES.txt', 'r') as cont_f:  #читаю названия файлов из похромосомных таблиц, склеиваю таблицы в одну
    for l in cont_f.readlines():
        l = '/home/ryabykh2018/all_to_all_data_4_12_22/9-backN2_DB/grid_nodup/mm10/ES/'+l.strip() # надо поменять путь, это не из базы, просто тестила
        i+=1
        if i == 1:
            contacts = pd.read_csv(l, sep='\t')
        else:
            contacts = contacts.append(pd.read_csv(l, sep='\t'), ignore_index=True)

chromosomes = pd.read_csv('mm10.chrom.sizes', sep='\t', header=None, names=['len'], index_col=0)
print(contacts)

'''bg calculating'''
contacts['bin'] = ((contacts.dna_start+contacts.dna_end+1000)/2000).round() # преобразую середину контакта к бину, из которого он
bg_contacts = contacts[(contacts.gene_type == 'protein_coding') & (contacts.rna_chr != contacts.dna_chr)]
contacts = contacts[['dna_chr', 'bin', 'count', 'gene_name_un']]
m = len(bg_contacts.gene_name_un.unique())
bg_contacts = bg_contacts[['dna_chr', 'bin', 'count']]
bg_contacts = bg_contacts.groupby(by=['dna_chr', 'bin']).count().reset_index()  # суммирую по бинам и рнк

for chr in bg_contacts.dna_chr.unique():
    bg_contacts.loc[(bg_contacts['dna_chr'] == chr), 'count'] = bg_contacts.loc[(bg_contacts['dna_chr'] == chr), 'count'] / bg_contacts.loc[(bg_contacts['dna_chr'] == chr), 'count'].sum()

for chr in chromosomes.index:

    unexist_bin = np.arange(chromosomes.len[chr]/1000) #добавляем ноль в каждый бин
    chr_l = [chr]*len(unexist_bin)
    cont_l = [0]*len(unexist_bin)
    zip_l = zip(chr_l, unexist_bin, cont_l)
    df = pd.DataFrame(zip_l, columns=bg_contacts.columns)

    bg_contacts = bg_contacts.append(df, ignore_index=True) 

bg_contacts.drop_duplicates(subset=['dna_chr', 'bin'], inplace=True) # убираю ноль, если есть не-ноль
bg_contacts.sort_values(by=['dna_chr', 'bin'], inplace=True)

for chr in bg_contacts.dna_chr.unique():
    bg_contacts.loc[(bg_contacts['dna_chr'] == chr), 'count'] = bg_contacts.loc[(bg_contacts['dna_chr'] == chr), 'count'].rolling(10,  min_periods=1,  center=True).sum() * 1/m * chromosomes['len'][chr]/10000 # усреднение по 10 бинов по процедуре из статьи 

bg_contacts = bg_contacts.reset_index()
bg_contacts.to_csv('grid_bg.tsv', sep='\t')

'''v-value calculating'''

contacts = contacts.groupby(by=['dna_chr', 'bin', 'gene_name_un']).count().reset_index()
result_df = pd.DataFrame(columns=['dna_chr',       'bin', 'gene_name_un',  'count'])

for rna in contacts.gene_name_un.unique():
    rna_contacts = contacts[contacts.gene_name_un == rna]
    print( rna_contacts)
    for chr in chromosomes.index:
        print(chr, rna)

        unexist_bin = np.arange(chromosomes.len[chr] / 1000)
        chr_l = [chr] * len(unexist_bin)
        cont_l1 = [rna] * len(unexist_bin)
        cont_l = [0] * len(unexist_bin)

        zip_l = zip(chr_l, unexist_bin, cont_l1, cont_l)
        df = pd.DataFrame(zip_l, columns=rna_contacts.columns)
        rna_contacts = rna_contacts.append(df, ignore_index=True)
        rna_contacts.drop_duplicates(subset=['dna_chr', 'bin'], inplace=True)
        rna_contacts.sort_values(by=['dna_chr', 'bin'], inplace=True)
        bg = bg_contacts.loc[(bg_contacts['dna_chr'] == chr), 'count'].reset_index() # добавление нулей в пропущенные бины (как выше)
        rna_contacts.loc[(rna_contacts['dna_chr'] == chr), 'count'] = rna_contacts.loc[(rna_contacts['dna_chr'] == chr), 'count']/ bg['count']/ \
                                                                      rna_contacts.loc[(rna_contacts['dna_chr'] == chr), 'count'].sum() * \
                                                                      chromosomes['len'][chr] / 1000 # подсчет обогащения


        rna_contacts['count'] = np.where(rna_contacts['count'] > 2, rna_contacts['count'], 0) # убиваем все обогащения, где меньше 2
    print(rna_contacts) # остатки от попытки найти место, где что-то идет не так
    for chr in chromosomes.index:
        rna_contacts.loc[(rna_contacts['dna_chr'] == chr), 'count'] = rna_contacts.loc[(rna_contacts['dna_chr'] == chr), 'count'].rolling(10, min_periods=1,  center=True).sum()/10 # усредняем по 10

    
        print(rna_contacts.loc[(rna_contacts['dna_chr'] == chr), 'count'])
    result_df = result_df.append(rna_contacts[rna_contacts['count']!= 0], ignore_index=True)
    result_df.to_csv('grid_method_result.tsv', sep='\t')

