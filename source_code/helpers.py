import requests
from Bio import Entrez
Entrez.email = 'stankovictamara.88@gmail.com'
import xmltodict
import re

import pandas as pd
import numpy as np
from Bio.Blast import NCBIWWW


def get_gen_bank_data(a):
    template = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={:s}&rettype=fasta&retmode=text'
    r_string = template.format(a['gene_id'])
    response=requests.get(r_string)
    whole_seq = ''.join(resp.text.split('\n')[1:])
    
    first = whole_seq[a['start'] - 10000:a['start']]
    second = whole_seq[a['end']:a['end'] + 10000]
    a['first']=first
    a['second']=second
    
    return a


def resolve_id(search_term):
    handle = Entrez.esearch(db="nuccore", term=search_term)
    t = xmltodict.parse(handle.read())
    return t['eSearchResult']['IdList']['Id']  

    
def fetch_by_id(gb_id):
    res_list = []
    handle = Entrez.efetch(db="nucleotide", id=gb_id, rettype="fasta_cds_na", retmode="text")
    f = handle.readlines()
    i = [line for line in f if line[0] == '>']
    
    for x in i:    
        primer=x.split(' ')[-2:]
        primer = [y.strip('[]') for y in primer]
        match = re.match(".*?(\d*)\.\.(\d*).*", primer[-1])
        prot_id = primer[0].split('=')[-1]
        beg = match.group(1)
        end = match.group(2)
        res_list.append({'protein_id': prot_id, 'beggining': beg, 'ending': end})
        df = pd.DataFrame(res_list)
        df.beggining = df.beggining.apply(np.int)
        df.ending = df.ending.apply(np.int)
    
    return df

def select_protein_IDs(df: pd.DataFrame, upstream, downstream):
    b1, e1 = upstream
    t1= df[(df.beggining >= b1) & (df.ending <= e1)]
    b2, e2 = downstream
    t2= df[(df.beggining >= b2) & (df.ending <= e2)]
    return (t1, t2)








