import xmltodict
import numpy as np
import pandas as pd




def processCrispr(crispr):
    return {'crisp_id': crispr['CRISPRId'], 
            'start': crispr['BeginningPosition'],
            'end': crispr['EndingPosition'],
             'dr': crispr['DR']['DRConsensus'],
            'hypothetical': crispr['Hypothetical']
           }

def generate_data_frame():

    raw_data = xmltodict.parse(open('data_from_crispr_finder/CRISPRdb.xml').read())
    raw_data = raw_data['CRISPRdb']['Taxons']['Taxon']

    rez = []
    for x in raw_data:
        name = x['ScientificName']
        seqs = x['Sequences']
        seqCount = int(seqs['SequenceCount'])
        if(seqCount == 1):
            crisprs = seqs['Sequence']['CRISPRs']
            cpCount = int(crisprs['CRISPRCount'])
            if cpCount == 1:
                d = processCrispr(crisprs['CRISPR'])
                d['name'] = name
                rez.append(d)
            elif cpCount > 1:
                for crispr in crisprs['CRISPR']:
                    d = processCrispr(crispr)
                    d['name'] = name
                    rez.append(d)
        else:
            for seq in seqs['Sequence']:
                crisprs = seq['CRISPRs']
                cpCount = int(crisprs['CRISPRCount'])
                if cpCount == 1:
                    d = processCrispr(crisprs['CRISPR'])
                    d['name'] = name
                    rez.append(d)
                elif cpCount > 1:
                    for crispr in crisprs['CRISPR']:
                        d = processCrispr(crispr) 
                        d['name'] = name
                        rez.append(d)
    df = pd.DataFrame(rez)
    df['end'] = df['end'].apply(np.int)
    df['start'] = df['start'].apply(np.int)
    df = df[['crisp_id', 'name', 'start', 'end', 'dr', 'hypothetical']]
    df['gene_id'] = df['crisp_id'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    result = df[df.hypothetical == 'No']
    result['blast_upstream'] = result.start.apply(lambda x: (x -10000, x))
    result['blast_downstream']=result.end.apply(lambda x: (x,x+10000))
    return  result

