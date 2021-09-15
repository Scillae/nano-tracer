import Reader, NanoConstructor
#诺诺帮忙改改

def DataReader_Nstar():
    dims_ls = [20,2,7]
    reader = Reader('data/6-arm-nanostar-starlike-412kInit.top',
                'data/6-arm-nanostar-starlike-412kInit.conf')
    strands_dic = reader.read_data()
    # Noella: Accepted_sci
    ns = NanoConstructor(strands_dic, dims_ls)
    return ns