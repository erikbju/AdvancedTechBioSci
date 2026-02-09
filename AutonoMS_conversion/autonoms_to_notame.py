import sys
import numpy as np
import pandas as pd

from pandas.io.formats import excel
excel.ExcelFormatter.header_style = None
excel.ExcelFormatter.index_style = None

# Function to convert the Mass spec data to a format that resembles MS-DIAL's output
def convert_report(input, metain, charge, align_id_start):
    if charge == 1:
        seq = "Positive"
        chrg = "POS"
    else:
        seq = "Negative"
        chrg = "NEG"

    # Removes all injections that are not QC, Sample, or Blank
    keep = metain[metain["Sequence"]==seq]
    keep = np.array(keep["Sample_Type"].isin(["QC", "SAMPLE", "BLANK"]))

    # Selects the features which were recorded during the given mode (Pos/Neg),
    # removes any duplicates, and selects the appropriate columns.
    # a = all features in the ymdb data base
    # b = all injections
    tmp = input[input["Precursor Charge"] == charge]
    a = tmp.drop_duplicates(subset= ['Molecule Name', 'Precursor Mz'], keep = 'first')
    a = a[["Molecule Name", "Collisional Cross Section", "Precursor Mz"]]
    b = tmp["Replicate Name"].drop_duplicates(keep = 'first').values

    # Assigns alignment_id for each feature
    alig_id = range(align_id_start, align_id_start + a.shape[0])

    # Creates a data frame where each row is a feature in the ymdb data base
    d1 = pd.DataFrame({"Alignment ID": alig_id, 
                       "Molecule Name": a["Molecule Name"],
                       "Cross Collisional Section": a["Collisional Cross Section"],
                       "Average Mz": a["Precursor Mz"],
                       "Column": ["IM"] * a.shape[0],
                       "Ion Mode": [chrg] * a.shape[0]})
    d1.reset_index(drop = True, inplace= True)    

    # Creates a matrix by iterating over the injections and extracting the 
    # area for each feature
    areas = np.empty((0,a.shape[0]))
    for i in b:
        c = tmp[tmp["Replicate Name"] == i]["Area"].values
        areas = np.append(areas, [c], axis = 0)

    # Converts NAs to zero, then transposes the matrix
    areas[np.isnan(areas)] = 0
    areas = np.transpose(areas)

    # Creates a dataframe based on the previous matrix, 
    # then removes injections that were not QC, sample, or blank.
    d2 = pd.DataFrame(areas, columns=b)
    d2 = d2.loc[:,keep]
    d2 = pd.concat([pd.DataFrame({"Alignment ID": alig_id}), d2], axis=1)
    d2.reset_index(drop = True, inplace = True)

    return d1, d2

# Function to generate a metadata data frame
def generate_meta(input):
    # Extracts data from the meta file
    inj_order = range(1,input.shape[0]+1)
    qc = input["Sample_Type"]
    Group = input["Description"]
    replicate = input["Replicate_Number"]

    # Creates data frame based on the data in the meta file
    out = pd.DataFrame({
        "Injection_order": inj_order,
        "QC": qc,
        "Group": Group,
        "Replicate": replicate
    })

    # Removes injections which are not QC, Sample, or Blank
    out = out[out.QC.isin(["QC", "SAMPLE", "BLANK"])]

    # Transpsoses the data frame
    outT = out.T
    return outT


def main():
    # Loads the files given by the terminal command
    meta = pd.read_excel(sys.argv[1],sheet_name="samples")
    pos = pd.read_table(sys.argv[2])
    neg = pd.read_table(sys.argv[3])
    
    # Generates the data frames, rows = aligned molecules and columns = injections
    posleft, posright = convert_report(pos, meta, 1, 0)
    negleft, negright = convert_report(neg, meta, -1, posleft.shape[0])

    # Combines positive and negative alignments
    outleft = pd.concat([posleft,negleft])

    # Adds the injections to the data frame
    out = outleft.merge(posright, on="Alignment ID", how="outer")
    out = out.merge(negright, on="Alignment ID", how="outer")

    # Removes duplicates
    out = out.drop_duplicates(subset = ['Average Mz', 'Cross Collisional Section'], keep='first')

    # Generates the meta data file
    metaout = generate_meta(meta)

    # The alignment data frame and meta data frame are combined into an excel file, 
    # to be used with notame
    fn = r'./output/comb_data.xlsx'
    with pd.ExcelWriter(fn) as writer:
        metaout.to_excel(writer, sheet_name='Sheet1', header=None, index=True, startcol=5, startrow=0)
        out.to_excel(writer, sheet_name='Sheet1', header=True, index=False, startcol=0,startrow=4)

if __name__ == "__main__":
    main()