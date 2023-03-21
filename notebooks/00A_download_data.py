#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
data_path = '/data2/hratch/FDA_collab/'
results_path = '/data2/hratch/FDA_collab/'


# In[ ]:


# directories
for dir_ in ['raw', 'raw/fda_data', 'raw/references', 'raw/references']:
    if not os.path.isdir(data_path + dir_):
        os.mkdir(results_path + dir_)
    
for dir_ in ['processed', 'interim', 'figures']:
    if not os.path.isdir(results_path + dir_):
        os.mkdir(results_path + dir_)



# In[ ]:


# lr pairs
if not os.path.isdir(data_path + 'raw/Ligand-Receptor-Pairs'):
    os.system('git clone https://github.com/LewisLabUCSD/Ligand-Receptor-Pairs.git ' + data_path + 'raw/')


# ### References

# Tabula Muris:

# In[31]:


gc1 = 'https://github.com/czbiohub/tabula-muris-vignettes.git'
os.system('git clone ' + gc1)

gc2 = 'https://github.com/czbiohub/tabula-muris.git'
os.system('git clone ' + gc2)

os.system('mv tabula-muris tabula-muris-vignettes/')
os.system('mv tabula-muris-vignettes ' + data_path + 'raw/references/')

hl = 'https://s3.amazonaws.com/czbiohub-tabula-muris/TM_facs_mat.rds'
cmd = 'wget -P ' + data_path + 'raw/references/tabula-muris-vignettes/data/ ' + hl
os.system(cmd)

hl = 'https://s3.amazonaws.com/czbiohub-tabula-muris/TM_droplet_mat.rds'
cmd = 'wget -P ' + data_path + 'raw/references/tabula-muris-vignettes/data/ ' + hl
os.system(cmd)


# fn = 'tabula-muris-vignettes/data/Data_Formatting.Rmd'
# cmd = 'Rscript -e "rmarkdown::render(' + "'" + fn + "'" + ')"'
# print(cmd)


# Brown DCs:

# In[5]:


# raw counts - downloaded directly from GSE137710 and copied into path
hls = ['https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4085nnn/GSM4085510/suppl/GSM4085510%5Fmouse%5Fspleen%5FTbx21CreRFPpos%5Fraw%5Fcounts%2Etsv%2Egz', 
      'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4085nnn/GSM4085511/suppl/GSM4085511%5Fmouse%5Fspleen%5FTbx21CreRFPneg%5Fraw%5Fcounts%2Etsv%2Egz']
# filenames: GSM4085510_mouse_spleen_Tbx21CreRFPpos_raw_counts.tsv.gz, GSM4085511_mouse_spleen_Tbx21CreRFPneg_raw_counts.tsv.gz
for hl in hls:
    cmd = 'wget -P ' + data_path + 'raw/references/ ' + hl
    os.system(cmd)
# metadata
hl = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137710/suppl/GSE137710%5Fmouse%5Fspleen%5Fcell%5Fmetadata%5F4464x9%2Etsv%2Egz'
cmd = 'wget -P ' + data_path + 'raw/references/ ' + hl
os.system(cmd)


# In[1]:


# raw counts - downloaded directly from GSE137710 and copied into path
hls = ['https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4085nnn/GSM4085510/suppl/GSM4085510%5Fmouse%5Fspleen%5FTbx21CreRFPpos%5Fraw%5Fcounts%2Etsv%2Egz', 
      'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4085nnn/GSM4085511/suppl/GSM4085511%5Fmouse%5Fspleen%5FTbx21CreRFPneg%5Fraw%5Fcounts%2Etsv%2Egz']
# filenames: GSM4085510_mouse_spleen_Tbx21CreRFPpos_raw_counts.tsv.gz, GSM4085511_mouse_spleen_Tbx21CreRFPneg_raw_counts.tsv.gz
for hl in hls:
    cmd = 'wget -P ' + data_path + 'raw/references/ ' + hl
    os.system(cmd)


# MCA:

# In[4]:


# # download MCA data
# hl = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108097/suppl/GSE108097%5FMCA%5FFigure2%5FBatchRemoved%5Fdge%2Etxt%2Etar%2Egz'
# cmd = 'wget -O ' + data_path + 'raw/references/ ' + hl
# os.system(cmd)

# # decompress
# fn = 'GSE108097_MCA_Figure2_BatchRemoved_dge.txt.tar.gz'
# cmd = 'tar xzf ' + data_path + 'raw/references/' + fn + ' -C ' + results_path + 'interim/'
# os.system(cmd)

# fn = 'MCA.zip'
# cmd = 'unzip ' + data_path + 'raw/references/' + fn + ' -d ' + results_path + 'interim/'
# os.system(cmd)


# Andreatta references:
# 
# https://www.nature.com/articles/s41467-021-23324-4#data-availability

# In[5]:


hls = ['https://ndownloader.figshare.com/files/23136746', 
      'https://ndownloader.figshare.com/files/23166794']

for hl in hls:
    cmd = 'wget -P ' + data_path + 'raw/references/ ' + hl
    os.system(cmd)

