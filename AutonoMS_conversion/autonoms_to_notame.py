import sys
import numpy as np
import pandas as pd

from pandas.io.formats import excel
excel.ExcelFormatter.header_style = None
excel.ExcelFormatter.index_style = None

def convert_report(input, metain, charge, align_id_start):
    if charge == 1:
        seq = "Positive"
    else:
        seq = "Negative"

    keep = metain[metain["Sequence"]==seq]
    keep = np.array(keep["Sample_Type"].isin(["QC", "SAMPLE", "BLANK"]))

    tmp = input[input["Precursor Charge"] == charge]
    a = tmp.drop_duplicates(subset= ['Molecule Name', 'Precursor Mz'], keep = 'first')
    a = a[["Molecule Name", "Collisional Cross Section", "Precursor Mz"]]
    b = tmp["Replicate Name"].drop_duplicates(keep = 'first').values

    chrg = "POS" if charge == 1 else "NEG"

    alig_id = range(align_id_start, align_id_start + a.shape[0])

    d1 = pd.DataFrame({"Alignment ID": alig_id, 
                       "Molecule Name": a["Molecule Name"],
                       "Cross Collisional Section": a["Collisional Cross Section"],
                       "Average Mz": a["Precursor Mz"],
                       "Column": ["IM"] * a.shape[0],
                       "Ion Mode": [chrg] * a.shape[0]})
    d1.reset_index(drop = True, inplace= True)    

    areas = np.empty((0,a.shape[0]))
    for i in b:
        c = tmp[tmp["Replicate Name"] == i]["Area"].values
        areas = np.append(areas, [c], axis = 0)

    areas[np.isnan(areas)] = 0
    areas = np.transpose(areas)

    d2 = pd.DataFrame(areas, columns=b)
    d2 = d2.loc[:,keep]
    d2 = pd.concat([pd.DataFrame({"Alignment ID": alig_id}), d2], axis=1)
    d2.reset_index(drop = True, inplace = True)

    return d1, d2

def generate_meta(input):
    inj_order = range(1,input.shape[0]+1)
    qc = input["Sample_Type"]
    Group = input["Description"]
    replicate = input["Replicate_Number"]

    out = pd.DataFrame({
        "Injection_order": inj_order,
        "QC": qc,
        "Group": Group,
        "Replicate": replicate
    })

    out = out[out.QC.isin(["QC", "SAMPLE"])]

    outT = out.T
    return outT


def main():
    meta = pd.read_excel(sys.argv[1],sheet_name="samples")
    pos = pd.read_table(sys.argv[2])
    neg = pd.read_table(sys.argv[3])
    

    posleft, posright = convert_report(pos, meta, 1, 0)
    negleft, negright = convert_report(neg, meta, -1, posleft.shape[0])

    outleft = pd.concat([posleft,negleft])

    out = outleft.merge(posright, on="Alignment ID", how="outer")
    out = out.merge(negright, on="Alignment ID", how="outer")

    out = out.drop_duplicates(subset = ['Average Mz', 'Cross Collisional Section'], keep='first')

    metaout = generate_meta(meta)

    fn = r'./output/comb_data.xlsx'
    with pd.ExcelWriter(fn) as writer:
        metaout.to_excel(writer, sheet_name='Sheet1', header=None, index=True, startcol=5, startrow=0)
        out.to_excel(writer, sheet_name='Sheet1', header=True, index=False, startcol=0,startrow=4)

if __name__ == "__main__":
    main()