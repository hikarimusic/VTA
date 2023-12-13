seq = '''
        1 aattccacaa ccttccacca aactctgcaa gatcccagag tgagaggcct gtatttccct
       61 gctggtggct ccagttcagg aacagtaaac cctgttctga ctactgcctc tcccttatcg
      121 tcaatcttct cgaggattgg ggaccctgcg ctgaacatgg agaacatcac atcaggattc
      181 ctaggacccc ttctcgtgtt acaggcgggg tttttcttgt tgacaagaat cctcacaata
      241 ccgcagagtc tagactcgtg gtggacttct ctcaattttc tagggggaac taccgtgtgt
      301 cttggccaaa attcgcagtc cccaacctcc aatcactcac caacctcttg tcctccaact
      361 tgtcctggtt atcgctggat gtgtctgcgg cgttttatca tcttcctctt catcctgctg
      421 ctatgcctca tcttcttgtt ggttcttctg gactatcaag gtatgttgcc cgtttgtcct
      481 ctaattccag gatcctcaac aaccagcacg ggaccatgcc ggacctgcat gactactgct
      541 caaggaacct ctatgtatcc ctcctgttgc tgtaccaaac cttcggacgg aaattgcacc
      601 tgtattccca tcccatcatc ctgggctttc ggaaaattcc tatgggagtg ggcctcagcc
      661 cgtttctcct ggctcagttt actagtgcca tttgttcagt ggttcgtagg gctttccccc
      721 actgtttggc tttcagttat atggatgatg tggtattggg ggccaagtct gtacagcatc
      781 ttgagtccct ttttaccgct gttaccaatt ttcttttgtc tttgggtata catttaaacc
      841 ctaacaaaac aaagagatgg ggttactctc taaattttat gggttatgtc attggatgtt
      901 atgggtcctt gccacaagaa cacatcatac aaaaaatcaa agaatgtttt agaaaacttc
      961 ctattaacag gcctattgat tggaaagtat gtcaacgaat tgtgggtctt ttgggttttg
     1021 ctgccccttt tacacaatgt ggttatcctg cgttgatgcc tttgtatgca tgtattcaat
     1081 ctaagcaggc tttcactttc tcgccaactt acaaggcctt tctgtgtaaa caatacctga
     1141 acctttaccc cgttgcccgg caacggccag gtctgtgcca agtgtttgct gacgcaaccc
     1201 ccactggctg gggcttggtc atgggccatc agcgcatgcg tggaaccttt tcggctcctc
     1261 tgccgatcca tactgcggaa ctcctagccg cttgttttgc tcgcagcagg tctggagcaa
     1321 acattatcgg gactgataac tctgttgtcc tatcccgcaa atatacatcg tttccatggc
     1381 tgctaggctg tgctgccaac tggatcctgc gcgggacgtc ctttgtttac gtcccgtcgg
     1441 cgctgaatcc tgcggacgac ccttctcggg gtcgcttggg actctctcgt ccccttctcc
     1501 gtctgccgtt ccgaccgacc acggggcgca cctctcttta cgcggactcc ccgtctgtgc
     1561 cttctcatct gccggaccgt gtgcacttcg cttcacctct gcacgtcgca tggagaccac
     1621 cgtgaacgcc caccaaatat tgcccaaggt cttacataag aggactcttg gactctcagc
     1681 aatgtcaacg accgaccttg aggcatactt caaagactgt ttgtttaaag actgggagga
     1741 gttgggggag gagattaggt taaaggtctt tgtactagga ggctgtaggc ataaattggt
     1801 ctgcgcacca gcaccatgca actttttcac ctctgcctaa tcatctcttg ttcatgtcct
     1861 actgttcaag cctccaagct gtgccttggg tggctttggg gcatggacat cgacccttat
     1921 aaagaatttg gagctactgt ggagttactc tcgtttttgc cttctgactt ctttccttca
     1981 gtacgagatc ttctagatac cgcctcagct ctgtatcggg aagccttaga gtctcctgag
     2041 cattgttcac ctcaccatac tgcactcagg caagcaattc tttgctgggg ggaactaatg
     2101 actctagcta cctgggtggg tgttaatttg gaagatccag cgtctagaga cctagtagtc
     2161 agttatgtca acactaatat gggcctaaag ttcaggcaac tcttgtggtt tcacatttct
     2221 tgtctcactt ttggaagaga aacagttata gagtatttgg tgtctttcgg agtgtggatt
     2281 cgcactcctc cagcttatag accaccaaat gcccctatcc tatcaacact tccggagact
     2341 actgttgtta gacgacgagg caggtcccct agaagaagaa ctccctcgcc tcgcagacga
     2401 aggtctcaat cgccgcgtcg cagaagatct caatctcggg aatctcaatg ttagtattcc
     2461 ttggactcat aaggtgggga actttactgg gctttattct tctactgtac ctgtctttaa
     2521 tcctcattgg aaaacaccat cttttcctaa tatacattta caccaagaca ttatcaaaaa
     2581 atgtgaacag tttgtaggcc cactcacagt taatgagaaa agaagattgc aattgattat
     2641 gcctgccagg ttttatccaa aggttaccaa atatttacca ttggataagg gtattaaacc
     2701 ttattatcca gaacatctag ttaatcatta cttccaaact agacactatt tacacactct
     2761 atggaaggcg ggtatattat ataagagaga aacaacacat agcgcctcat tttgtgggtc
     2821 accatattct tgggaacaag atctacagca tggggcagaa tctttccacc agcaatcctc
     2881 tgggattctt tcccgaccac cagttggatc cagccttcag agcaaacacc gcaaatccag
     2941 attgggactt caatcccaac aaggacacct ggccagacgc caacaaggta ggagctggag
     3001 cattcgggct gggtttcacc ccaccgcacg gaggcctttt ggggtggagc cctcaggctc
     3061 agggcatact acaaactttg ccagcaaatc cgcctcctgc ctccaccaat cgccagtcag
     3121 gaaggcagcc taccccgctg tctccacctt tgagaaacac tcatcctcag gccatgcagt
     3181 gg
'''

seq = seq.split()
seq = [s for s in seq if not s.isnumeric()]
seq = ''.join(map(str, seq))
for i in range(len(seq)//200):
    print(i)
    print("map", seq[200*i:200*i+100], seq[200*i+100:200*i+200], '\n')
    print("")
# print(seq)