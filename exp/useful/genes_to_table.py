import math

def create_latex_table(data, cols=6):
    # Split gene names from the comma-separated data
    genes = [line.split(',')[0] for line in data.strip().split('\n')]
    
    # Calculate rows needed
    total_genes = len(genes)
    rows = math.ceil(total_genes / cols)
    
    # Start LaTeX table
    latex = [
        "\\begin{table}[htbp]",
        "\\centering",
        "\\caption{Caption}",
        "\\small",
        "\\begin{tabular}{" + "l" * cols + "}",
        "\\hline"
    ]
    
    # Add headers
    headers = [f"\\textbf{{Column {i+1}}}" for i in range(cols)]
    latex.append(" & ".join(headers) + " \\\\")
    latex.append("\\hline")
    
    # Fill table by columns
    for row in range(rows):
        row_genes = []
        for col in range(cols):
            idx = row + col * rows
            if idx < total_genes:
                row_genes.append(genes[idx])
            else:
                row_genes.append("")
        latex.append(" & ".join(row_genes) + " \\\\")
    
    # Close table
    latex.extend([
        "\\hline",
        "\\end{tabular}",
        "\\label{tab:label}",
        "\\end{table}"
    ])
    
    return "\n".join(latex)

# Sample usage (replace with actual data reading)
data = """HEG1,1.149956745258744,6.425230752464496e-36,1.0213546804117563e-31
ADGRA2,1.796381635613316,7.294227310951484e-34,5.7974518667442396e-30
ANTXR1,2.0633284849173714,4.331521543707271e-30,2.250086968679463e-26
PCDHGC3,1.8327916523099663,5.662020555308161e-30,2.250086968679463e-26
F2R,1.3800401270662592,4.2790408627682035e-29,1.1336605592427229e-25
ETS1,1.0619200849816621,1.9875503875679031e-28,4.513442994397056e-25
KIRREL1,1.7136044020081296,2.5461886807186346e-27,5.059276908587927e-24
KCTD12,1.40893286886319,1.2141367388546451e-26,2.144435288981493e-23
CD93,1.073690930530993,2.099090232952401e-26,3.3367138343011365e-23
GEM,2.1315921460288094,3.165336361790869e-26,4.5741988006388776e-23
FBN1,1.6501259512509638,9.174826949201173e-26,1.215358743204182e-22
VASH1,1.0103331916245073,1.0867021147321045e-25,1.3287859089062719e-22
GUCY1A1,1.8052372106470038,2.0192190342278033e-25,2.215106347086905e-22
DCN,2.3462358521776254,2.0902488177090825e-25,2.215106347086905e-22
IQGAP1,1.2535787453622271,3.5326997866413094e-25,3.509737238028141e-22
AP001189.3,1.8938616117593674,4.156621938536126e-25,3.8866860197041326e-22
HSPG2,1.4047557927827463,1.0796069117501014e-24,9.486494318164445e-22
FNDC1,2.7024616862480535,1.1338914949995248e-24,9.486494318164445e-22
SPRED1,1.0681203259081136,1.6723956804823082e-24,1.3292200868473386e-21
THBS2,2.799704019613927,1.9437635512263274e-24,1.4713364481092239e-21
GUCY1B1,1.2187522827921486,5.031283774828541e-24,3.477273342811935e-21
LUM,2.491673313689831,1.0704372982797308e-23,6.80626851738184e-21
CTTNBP2NL,1.1055099475207237,1.4016829282997056e-23,8.569673780096969e-21
MSRB3,1.2837723622342623,1.891434044616175e-23,1.1135642804895823e-20
SYT11,1.0632453196591736,2.983987366030541e-23,1.6940522560864815e-20
HIF1A,1.2429243218212094,3.225156405276342e-23,1.767830559250784e-20
ZNF532,1.1816744800860304,9.069044366170644e-23,4.4574683990678044e-20
SLC4A7,1.0201338093950993,1.031627192370412e-22,4.6904045601982554e-20
ASPN,1.489688206645716,1.3204892366040326e-22,5.673107271637216e-20
ADAMTS1,1.523208670826968,1.8146434339341624e-22,7.475982984377563e-20
COL12A1,1.6717518632897983,2.155806131841866e-22,8.358218115062999e-20
ADCY3,1.3653952198682293,2.3104461795814093e-22,8.7444886834824e-20
RAB31,1.287930869704184,2.9705521870274714e-22,1.0265195122823628e-19
TEK,1.4553675339817134,3.090084185347887e-22,1.0451059193678725e-19
PHLDB1,1.5049171703419721,4.141135210253326e-22,1.37140594379556e-19
CRISPLD2,1.896556325612367,8.841005540160856e-22,2.755620079733274e-19
PDE1A,1.7035040680241211,1.222916088324933e-21,3.6678253094364406e-19
COL4A1,1.0856716343251862,4.500622316232022e-21,1.2330788129442988e-18
FAP,2.4519335802153392,5.182723537956464e-21,1.3505667763828845e-18
DSE,1.456822864028112,5.735534153425015e-21,1.4471754111562545e-18
PODN,2.149547557350786,5.9617995900792286e-21,1.4807619731859285e-18
EPHA3,2.364204579003555,7.819062499347083e-21,1.8832093559033522e-18
PEAR1,1.0061076041073718,8.143709699927905e-21,1.9321255132843878e-18
PXDN,1.3717222142273853,1.1412469947936349e-20,2.6678326807705322e-18
THBD,1.3477557243980265,1.1769623421071278e-20,2.7114483174108557e-18
CALHM2,1.0090600282265323,1.204758093502925e-20,2.7358335220460707e-18
COL5A1,2.015023749399842,1.2322914317657862e-20,2.7589443097674564e-18
ARHGAP23,1.129409632508855,2.023513013071028e-20,4.406268884353022e-18
COL8A1,1.7018769200000063,2.45520778898658e-20,5.203731068497424e-18
HEPH,1.533935967449107,2.5075858250567235e-20,5.2448137204081155e-18
PTPRE,1.4847039566724498,2.796096499301545e-20,5.6983012760124816e-18
LTBP4,1.2895898145133151,3.321809838927224e-20,6.683985974631285e-18
ARHGAP31,1.0783461883834602,3.864331384766678e-20,7.61076438926134e-18
EDNRA,1.3014160406862407,3.878157495786163e-20,7.61076438926134e-18
ITGA9,1.3993742770532558,4.316965859528602e-20,8.368596256471543e-18
ZEB2,1.4349224601585295,7.481786061365055e-20,1.399182014487752e-17
PLSCR1,1.0176075458887008,8.086427347333733e-20,1.4946726641071747e-17
FOXF1,1.6023637052400046,8.376652138804885e-20,1.5131279818004824e-17
AXL,1.270222165028825,9.533109369814651e-20,1.665256115852458e-17
COL4A2,1.1189874614848276,1.0665955594351394e-19,1.823075592772148e-17
SVEP1,2.124859594652141,1.152731353244933e-19,1.9493422969341973e-17
ITGA11,1.4331774274760305,2.0249815322238954e-19,3.2846026975745957e-17
JAM3,1.490139828179416,3.0529869684372235e-19,4.8530280850278104e-17
KLF7,1.0916331436275886,3.5927960572786016e-19,5.654562982821846e-17
FILIP1L,1.1728139221192362,7.358431381712292e-19,1.0537804076008882e-16
SPARC,1.0147170772021221,8.133587246123641e-19,1.1543884184319768e-16
NEXN,1.488815568217822,8.448242832918339e-19,1.1858444792967608e-16
NOTCH3,1.2432829000794163,9.548898193771248e-19,1.3085283248981704e-16
FRMD6,1.0920019264208949,1.0658897642259e-18,1.435879973909738e-16
LAMA2,1.894714640023317,1.165227373252255e-18,1.556508766825029e-16
TNXB,2.0520926755979225,1.7837616501108468e-18,2.232651589776537e-16
DOCK10,1.223465766017081,2.161689650696024e-18,2.623070128814046e-16
ECM1,1.7962044716538654,2.189145422262194e-18,2.63626179032423e-16
IFI16,1.231385122369764,2.3288570915474565e-18,2.7421860983139533e-16
FPR3,1.264748331016693,4.609129683978938e-18,4.98413098343736e-16
LTBP2,1.7257027946488066,5.716633104413865e-18,5.831532202592026e-16
OSMR,1.3256608266021026,6.177544003494595e-18,6.21507844807279e-16
MAML2,1.1166493610879082,6.537297548549405e-18,6.516486343178728e-16
ZSWIM4,1.422210066960845,6.961635936800128e-18,6.873426388284151e-16
DCHS1,1.1126019378728684,7.241677797408633e-18,7.105784584420225e-16
MYOF,1.903000356883081,8.549586696180352e-18,8.286843300151394e-16
SNAI1,1.357510230191282,1.0917816988244894e-17,1.0269208215688806e-15
LZTS1,1.080723235544498,1.4250686821080294e-17,1.2937013822660943e-15
EMP1,1.4325406697841747,1.432382003515555e-17,1.2937013822660943e-15
ADAMTS2,1.7308787982379028,1.7352226796247872e-17,1.5409552913584143e-15
COL3A1,2.2697840230115087,2.1122226570853542e-17,1.8448291954411422e-15
SSC5D,1.960234417304748,2.9057387945356675e-17,2.4568948871244133e-15
COL1A2,2.0822988110642338,3.387015152405307e-17,2.7896369358878115e-15
MEDAG,1.9839529903627253,5.829727531584394e-17,4.564992553796331e-15
ARHGEF17,1.0468260316857996,6.051233942516562e-17,4.715216409325651e-15
CCN1,1.3154348289928466,8.457144687542862e-17,6.341262827980251e-15
NCOA7,1.0860154750721434,9.197573576776781e-17,6.831992036282417e-15
SLIT2,2.1735177960694636,1.0071136982848987e-16,7.377455920708179e-15
LBH,1.0480609436991357,1.0199517943479239e-16,7.431185336981756e-15
CLEC7A,1.2646138473498654,1.0237981811770285e-16,7.431185336981756e-15
DACT3,1.3778579558173782,1.0506566320920563e-16,7.591471738061512e-15
CYBRD1,1.2536634086448692,1.0646046993388978e-16,7.657446289905483e-15
AMPD3,1.0584533398052178,1.147541077737979e-16,8.216807644920232e-15
BGN,1.2708684312854395,1.3984524033266788e-16,9.749912018982844e-15
PCDH18,1.418078856280066,1.6704441452810577e-16,1.1347598347601578e-14
PLXDC2,1.3218085208647263,1.9477434761261481e-16,1.2793938139050104e-14
PLAT,1.4031695447933046,1.9727761914528001e-16,1.2905041291906877e-14
ICAM1,1.2349206601757674,1.9853755273321478e-16,1.2934233353472059e-14
PRR16,1.2943649327865303,2.0270412243008761e-16,1.3151774408770094e-14
CPXM2,1.466930526820287,2.1270697605806064e-16,1.3689028710198106e-14
COL1A1,2.100853892344803,2.5704997454690804e-16,1.6086875572431692e-14
NGFR,1.6215891392456068,2.700770566578081e-16,1.670484394020435e-14
ITGBL1,2.108586448480734,3.03260553306804e-16,1.8469845805996004e-14
ITGAV,1.0346365761636713,3.0503458396837056e-16,1.8506983766264195e-14
RCAN2,1.7426204758527597,3.4348718027184837e-16,2.0526587284215418e-14
ITGA4,1.2021323319839596,3.958546618213053e-16,2.330557668263507e-14
AC090559.1,1.0297022074911406,4.0977887540930554e-16,2.4036328426222585e-14
DACT1,1.3931291888856934,4.641676654474179e-16,2.6636856353617887e-14
CCDC80,2.2960650467726964,4.720617984620298e-16,2.6992425713497933e-14
MOXD1,2.1665152303384363,4.909936936002334e-16,2.767672253003302e-14
CXCR4,1.0125556890108958,6.843396288741088e-16,3.7028565317176156e-14
EGR2,1.9272245046805803,7.396041193403805e-16,3.9718740138630704e-14
GAS1,2.163752975249805,7.561409767090757e-16,4.033428512002506e-14
SLC9A9,1.1004708280261828,7.686245160795465e-16,4.086306122943302e-14
TRERF1,1.001271112829237,7.823582749829748e-16,4.136611462227838e-14
ABI3BP,1.6175308771381511,7.832914255979991e-16,4.136611462227838e-14
SULF1,2.0210230606987647,8.000595117992222e-16,4.197275907445689e-14
EFEMP1,2.0432495805200395,1.0448214423483062e-15,5.340347796645876e-14
GPR183,1.2951796028895424,1.7094444095280137e-15,8.386829732672008e-14
PTGIR,1.199551648684725,1.8747936135798192e-15,9.08589002483683e-14
POSTN,2.3820625862569043,2.0491876464730637e-15,9.811411695281873e-14
YWHAZP5,1.0333109253424941,2.303622404660703e-15,1.0833840752806668e-13
C7,2.118241720799431,2.787604724466202e-15,1.2918881836768148e-13
FLNA,1.4853945374598054,2.9210015654325214e-15,1.3304367015505834e-13
OSBPL5,1.0507101742026272,2.955713578123319e-15,1.3424006582242364e-13
EHD3,1.1798162623061912,3.015300392937755e-15,1.3655616822261696e-13
FMOD,2.0194391850092592,3.254261642729512e-15,1.4612921771985402e-13
LXN,1.4872754322156783,3.3802164477384077e-15,1.5135752296690063e-13
CXCL12,1.5552207950958634,4.287399266750597e-15,1.8417571272832245e-13
CD200,1.0474173091052543,5.4603954437891815e-15,2.296255184509863e-13
PDGFRA,2.0924381265953067,8.218534110022337e-15,3.2660454553228764e-13
COL14A1,1.8352520890231292,9.9462784607604e-15,3.8189865316967946e-13
ID4,1.4753542659232042,9.991566044611044e-15,3.827130936027401e-13
EMILIN1,1.7097028219505137,1.2753908144033352e-14,4.73682532377463e-13
HEYL,1.1851963555215361,1.4462340573062158e-14,5.284904959756231e-13
PTGIS,2.3557915330930355,1.542927224486966e-14,5.573231512008049e-13
AC015922.2,1.567911537264084,2.6952040300108244e-14,9.000622533834467e-13
DDR2,1.2900392438761943,3.341504255864799e-14,1.0929331615478774e-12
MECOM,1.3237261244908822,3.46445924241068e-14,1.121609859824036e-12
MYH11,1.4586843878208124,9.74130721286923e-14,2.780032665274134e-12
TCIM,1.2719521261464009,1.0754872075001933e-13,3.0419830338831092e-12
FPR1,1.6959477333601345,1.1066893368058464e-13,3.1191371804726483e-12
MMRN1,1.3078120202008512,1.4738201075281398e-13,4.018498186838304e-12
NID2,1.3282816964282838,1.513747405662177e-13,4.0992382896773365e-12
OSM,1.5327299936586658,1.720417219208483e-13,4.619552722388183e-12
LAMB1,1.0246129966281698,1.7245979878336034e-13,4.622969580877396e-12
ARSJ,1.271398173094125,1.9149927402452957e-13,5.090422173735655e-12
SERPINB9,1.0004206445681814,1.9324207079101013e-13,5.1196265954898285e-12
CDR2L,1.6111967162814067,1.9502477504517746e-13,5.158259274739003e-12
UNC5B,1.1664488722402526,2.025664361637273e-13,5.348830679831577e-12
MMP28,1.423201606713407,2.127314580022696e-13,5.598641152986883e-12
IL1B,1.5861313146421914,3.0168747623365164e-13,7.636344143646697e-12
CSF2RB,1.1425251511941357,3.3300916222798373e-13,8.336241957127606e-12
GFPT2,2.244843414764091,3.528745554293211e-13,8.791996760351863e-12
THBS1,1.5051558480722895,3.711209923055965e-13,9.117989634760066e-12
AKT3,1.3649482806203037,3.8279939878116926e-13,9.347126333372454e-12
FAM171B,1.3318241017620995,4.3837189299765905e-13,1.0510346321403903e-11
CHST11,1.1566827831685171,4.887837133175231e-13,1.1579293452899177e-11
STX11,1.0576766428551012,6.043689367903411e-13,1.3923258868433715e-11
DOCK8,1.0442224062275824,9.798758114051048e-13,2.1484283997373166e-11
ENDOD1,1.340453879808756,9.907009393035936e-13,2.166187363297101e-11
LRRC32,1.2274424377091342,1.0472112408197089e-12,2.2772188623898895e-11
TMEM119,2.074354450757444,1.0641024486947618e-12,2.3076360879197728e-11
ZFP36,1.065565303903002,1.0759966515430747e-12,2.3270806493780566e-11
COL6A3,1.9034937523453355,1.1675165075435902e-12,2.5045671260341306e-11
CALCRL,1.0154221621053208,1.3166335202520625e-12,2.775756822006205e-11
STING1,1.1897307761348697,1.4453676197494337e-12,3.001696304821773e-11
F3,1.6909265014322137,1.5826750361461814e-12,3.254618677177193e-11
RGS1,1.4564378147181198,1.5915247676057972e-12,3.2685888508865315e-11
SLC25A24,1.210070572347063,1.6987759400101636e-12,3.4487538112901096e-11
SLC6A6,1.948570734147107,1.8306553893094547e-12,3.660389694146301e-11
SOCS3,1.912059695754368,1.9277783154030404e-12,3.825713370992101e-11
PFKFB3,1.8061833204867068,1.9520390388287745e-12,3.8594045474156963e-11
MNDA,1.16874492605005,2.608635660024686e-12,4.9483141350539864e-11
HMCN1,1.248596960280216,2.8052380613283685e-12,5.2771673636539346e-11
MMP19,1.0234332206576875,3.0676891802518195e-12,5.6888012333875834e-11
COL10A1,3.532274608771236,3.1879421948936294e-12,5.872019597917629e-11
RND3,1.2928060154285754,3.728801343253397e-12,6.775361857979796e-11
DOCK2,1.0673364408604273,3.9503951500002975e-12,7.13585014822781e-11
CDH11,2.2349874743680247,3.96816163682157e-12,7.151689045228534e-11
B3GNT9,1.2432685537182349,4.007085588369363e-12,7.213661666219636e-11
PMP22,1.3482229907936374,4.189972386856752e-12,7.479830725008439e-11
HPSE,1.1298344304601438,4.443162630203882e-12,7.85634184312802e-11
LIFR,1.6750479781723513,6.004164767822243e-12,1.0273649424036854e-10
CCN3,1.3986895764544038,6.1168791841993385e-12,1.0423700947254552e-10
COL15A1,1.1573954236543882,6.40456940129638e-12,1.0865211867983699e-10
ALOX5,1.873352902411173,6.506547821296389e-12,1.1014705449129649e-10
UBASH3B,1.4054395783956664,7.010675233978157e-12,1.1743065702773106e-10
DPYSL3,1.5071354899621205,8.82654505257699e-12,1.43904369390527e-10
SMOC2,1.3176297199609142,8.851884516490911e-12,1.44169627330061e-10
PRELP,2.2526986736152645,8.991342927971112e-12,1.459922238846055e-10
CHSY1,1.0125481073491955,9.643478474747848e-12,1.554693040918781e-10
COL16A1,1.708605045465017,1.039815308891751e-11,1.652177021336877e-10
GAL3ST4,1.530383422888001,1.1143288590544411e-11,1.7503331564752367e-10
VCAN,2.671474485474432,1.2039977911107315e-11,1.8818828797931352e-10
FBLN5,1.1889151553099002,1.2991363966995533e-11,2.010815205641295e-10
ITPR3,2.089162740263632,1.475577464796548e-11,2.2445721895125287e-10
SAMSN1,1.0769468981593513,1.5071683294202177e-11,2.2797350591494335e-10
COL6A2,1.5254486939052514,1.6196423298826444e-11,2.428852309039105e-10
ITGA3,2.3565994000218238,1.666705823971627e-11,2.485361705239492e-10
SAMD9,1.0227616331844012,1.7290693450673645e-11,2.5663199168245403e-10
PRSS23,1.1558153444334343,1.8589967706465413e-11,2.7260712791695034e-10
SAMD9L,1.0747129265943658,1.976059040594274e-11,2.8765049916929104e-10
MFAP4,2.169055751057289,1.9862729022461888e-11,2.888727726816598e-10
GPX8,1.5052805839513799,2.296794464377419e-11,3.283259424976929e-10
SIGLEC7,1.194775686009473,2.350811075204044e-11,3.3454335587684404e-10
SEMA3C,2.394251886956661,2.4802892140311946e-11,3.510835026379329e-10
KLF4,1.1766555946004256,2.6485057914991654e-11,3.72572106740449e-10
LIMCH1,1.384971993797404,2.6656449512677548e-11,3.7432060199074414e-10
SPON1,1.9015380598215812,2.9642990197076693e-11,4.1080318960754536e-10
TLR2,1.083356262421815,3.038345146636973e-11,4.1961367898298285e-10
SSPN,1.8485050407829444,3.2560807160852e-11,4.469234870472461e-10
MB21D2,1.0023842661034525,3.8251220416939506e-11,5.166027185621669e-10
TLR1,1.0501767229918901,4.0106048244831455e-11,5.370899266216013e-10
PTGS1,1.306144302748763,4.4206944843449926e-11,5.831648093207302e-10
MS4A6A,1.001364519313236,4.872355479112601e-11,6.369322590129434e-10
TMEM159,1.480483083969243,5.560817136001843e-11,7.198269478329421e-10
MGP,1.2907469878764413,6.481547052694759e-11,8.222719229819304e-10
TAGAP,1.1256399507092467,6.543691220440726e-11,8.288327939452254e-10
CD163,1.2090427571882003,6.810246102468273e-11,8.591720003558386e-10
FAR2,1.0509689051063373,7.508310065168668e-11,9.338974710165975e-10
NFASC,1.037228995409921,9.784917548748814e-11,1.1821614982113525e-09
ELF4,1.1240857843941852,1.1476034412533484e-10,1.3644206658312062e-09
ATP1A1,1.1147396134169834,1.1630835590701557e-10,1.3807599891694696e-09
ABCC1,1.1156657894375916,1.2394807951906628e-10,1.4573067100851165e-09
GLIS3,1.4546623661585405,1.2604754668734773e-10,1.4808956409032367e-09
GGT5,1.44546650283517,1.311155720706798e-10,1.5347666668891945e-09
TGFB3,1.4612864878393834,1.3127888294518395e-10,1.5355475520946608e-09
FIBIN,1.9904403395860226,1.3783831638109518e-10,1.6063620800541708e-09
PAPLN,1.8324847789150625,1.6665364559467102e-10,1.8972646986675998e-09
IL7R,1.6851236592526855,2.1971326668148412e-10,2.4287636211188257e-09
ADAMTSL2,1.1287111887256693,2.200833125405247e-10,2.4299010366905696e-09
SULF2,1.104062715365189,2.455471183097293e-10,2.6734362963366146e-09
SYNPO2,1.500409274517014,2.831732456456087e-10,3.0169717914092466e-09
CHSY3,1.0270436197869366,2.8539061767908434e-10,3.0340915928699094e-09
P2RY13,1.1645776459831152,3.303768200350787e-10,3.4482402700443934e-09
FYB1,1.1174279101857514,3.345792383811196e-10,3.480675113420339e-09
CD69,1.1189136757266631,4.251835185816693e-10,4.324195272792205e-09
RAB11FIP1,1.251374863652861,4.6313542135586543e-10,4.680229280275167e-09
PTGER4,1.0519575429238615,4.724506010316395e-10,4.753211869619584e-09
LPAR5,1.1628740853015607,5.963553958648014e-10,5.815745627402995e-09
GABBR1,1.676533476140969,7.04466580812153e-10,6.669565675157823e-09
SIGLEC10,1.243231309314311,7.682185973178509e-10,7.195994592200681e-09
FZD1,1.5484699840465823,7.925242881366466e-10,7.38450532486526e-09
KIF3C,1.7097922465148896,8.380705518304469e-10,7.75886400226953e-09
IER3,1.2168397326789966,1.0538099412501735e-09,9.476406079891987e-09
PELI2,1.4385103978546294,1.0666373805004365e-09,9.563038804531833e-09
MARVELD1,1.0617205134513779,1.14127989430171e-09,1.016916210752241e-08
FOS,1.6718002135175447,1.3017312910213985e-09,1.1432221327113896e-08
LAMC3,1.1563547913708578,1.3643602460556427e-09,1.1896802233297036e-08
JAG1,1.1002147167578145,1.409364268214878e-09,1.2248908916098252e-08
ANXA1,1.236904249669187,1.6445034102131745e-09,1.4031683418544616e-08
INHBA,1.1100966515993984,1.7079600723472677e-09,1.4503062665615474e-08
SHISA3,1.7354342103730114,2.1641653325408388e-09,1.7824648769984027e-08
ALOX5AP,1.2006687307900985,2.367346823531923e-09,1.9367650595400644e-08
CXCL8,1.8790045493154406,2.396204491064486e-09,1.957351828877753e-08
COL4A4,1.1031078151272677,2.41415558152809e-09,1.9659537461050472e-08
FGFR1,1.4898014171981804,2.4197165850340083e-09,1.9689210773919167e-08
CELF2,1.3942461801385977,2.48190211194459e-09,2.0108214052737617e-08
CCL2,1.7915629001786626,2.6497878067059122e-09,2.129475580151526e-08
SLIT3,1.4199244670037163,2.7600251942347593e-09,2.2013728292802674e-08
FUT4,1.3025736042208318,2.8750024500300594e-09,2.2861950448062943e-08
CTBP2,1.1334109058238457,2.936982794087968e-09,2.3250138692640606e-08
F2RL1,1.0510594663289594,3.3988225255994293e-09,2.634211743877549e-08
SPART,1.0625281484984312,3.613955955781622e-09,2.7806119977301386e-08
SLC12A2,1.0639492708439093,3.851295464319704e-09,2.94327849523202e-08
OGN,2.1284388613999785,3.8655310509202565e-09,2.9513199608755232e-08
PTPRC,1.227282473143242,4.4668124755964714e-09,3.353385583524753e-08
SPIRE1,1.0982554116326382,4.5205543922348956e-09,3.3879647627989585e-08
MXRA5,2.4856072206967816,5.3385198057754055e-09,3.930574841714027e-08
NCEH1,1.7528392235130439,5.411964508435388e-09,3.981801180819471e-08
ARMCX2,1.5692921897431338,5.54457063423207e-09,4.069090249388412e-08
EPB41L3,1.2240542738114588,5.579922771912981e-09,4.0912570287052e-08
HTRA3,1.3411038618993016,5.779828683179817e-09,4.2280790035815166e-08
SMIM10,1.127172085804115,5.902179895755164e-09,4.305693052910697e-08
GPRC5B,1.0691948084761986,5.908865578848314e-09,4.308592992723523e-08
LHFPL6,1.390081135764352,6.052557964810095e-09,4.3992437772574885e-08
TUBB6,1.092677033737958,6.686683243817818e-09,4.803050919282785e-08
RIPK3,1.0292253505110611,7.066318452977435e-09,5.050638405059771e-08
PTAFR,1.0220800729395705,8.334097358766704e-09,5.85153761550157e-08
C11orf96,1.4820011524193109,8.523870154272325e-09,5.971592770926086e-08
VSIG4,1.3826871654276165,8.836534972491211e-09,6.168887128797553e-08
ADCY5,1.2654456778462562,9.27201377034346e-09,6.424931599537037e-08
EGR1,1.1833120941594752,1.0379193566306196e-08,7.081286668289943e-08
NPY1R,1.6677583825853568,1.1850915853580319e-08,7.989065242091296e-08
LTBP1,1.0049938776717355,1.3845420118826493e-08,9.166463898744936e-08
FSTL1,1.3553826790881898,1.4060385088816576e-08,9.293300680741302e-08
SEL1L3,1.7035717661607168,1.627631241669147e-08,1.0534538362203892e-07
CYBB,1.0921647451562178,1.742958933380652e-08,1.1207959225331247e-07
PDPN,2.273041783953413,1.843027608918272e-08,1.1779962553825836e-07
COLEC10,2.2142485424231286,2.003898242354425e-08,1.272631500617896e-07
EDN1,1.1710475065921422,2.1511551578932723e-08,1.3556960811211976e-07
LPAR1,2.0705307051909605,2.1517496305163448e-08,1.3556960811211976e-07
APOBEC3C,1.0697520690240019,2.2988068514904088e-08,1.4403560784900093e-07
GAS7,1.1021491041032563,2.532372140627705e-08,1.574289696809464e-07
SFRP4,2.0358502136502086,2.6633226808240595e-08,1.6466813432275088e-07
ISLR,2.3674285364833727,3.1328989961190266e-08,1.9080675265252125e-07
TPBG,1.7740766193709987,3.7270703980194406e-08,2.2475535298526945e-07
KCNE4,1.9770128032699923,4.1321956820572184e-08,2.465667513587896e-07
TIMP2,1.4462714854506276,4.1725749793449784e-08,2.488827462351511e-07
PDE5A,1.302852495869217,4.655942749982683e-08,2.7388915687408404e-07
COMP,2.2561422933878355,4.8062573643903876e-08,2.821280172243338e-07
BCO2,1.4788156556512286,4.9568161115322245e-08,2.904295942090536e-07
IL2RA,1.2393504481657553,5.7077181716492045e-08,3.300468827083876e-07
CXCL3,1.717241471207835,6.095132411885509e-08,3.5142627790834983e-07
PGM2L1,1.1230269422502257,6.485294802653831e-08,3.720326459147791e-07
CARD11,1.1134485136352208,6.681772535016904e-08,3.819254089055329e-07
PLD4,1.2229208630537234,7.121166415726621e-08,4.054371824655815e-07
RRAD,1.7423000937020896,7.174957929002356e-08,4.082073415870489e-07
CSGALNACT1,1.3656606685440154,8.508309755064088e-08,4.767292628357375e-07
CCDC9B,1.098060857699175,1.0466871191979402e-07,5.781146089913293e-07
ABR,1.4588037617559042,1.0701297476380261e-07,5.896285084386157e-07
DZIP1,1.158976245666561,1.1228658654963678e-07,6.148493213203673e-07
COL6A1,1.104982662468299,1.1311424249921904e-07,6.187419128587701e-07
GDF2,3.551343207522356,1.1664546623242051e-07,6.363062221106919e-07
MYO10,1.5806622371254877,1.3607532805849562e-07,7.273608544290371e-07
SERPINE1,1.3480161953286904,1.4603039197595235e-07,7.750581338396456e-07
INMT,1.4596097246394448,1.6801507887572394e-07,8.802793980911365e-07
GLT8D2,1.6327220170258414,1.727036689583022e-07,9.024646685605431e-07
ADGRG1,1.1091463096166696,1.8212835701808121e-07,9.464244403920951e-07
MMP7,2.57304331018206,1.8719485670387595e-07,9.689512999559792e-07
SFRP1,2.7860302745917185,1.8771165697779766e-07,9.713100583720937e-07
PTH1R,1.4830616329164874,1.9384491973242216e-07,1.001091242386804e-06
NR4A3,1.4462781631843016,2.3114444962373138e-07,1.1735139480098479e-06
CHST1,1.1583597065462226,2.378437857036478e-07,1.2041129965845177e-06
LDLRAD3,1.0698458514588853,2.4392062281133444e-07,1.2301276079343186e-06
GALNT7,1.4854549889064002,2.545166000768914e-07,1.2782925354888676e-06
SCARA5,2.937770159126445,2.6039189858793697e-07,1.3020414029423862e-06
PMEPA1,2.027652148325915,2.6241232445182907e-07,1.3109070740057434e-06
TGFA,1.5257402152660873,2.7704001246128153e-07,1.3753366764786169e-06
MAP1B,1.287519480028709,2.9901865913879423e-07,1.475720048955901e-06
RAI2,1.0530481380619712,3.010474593267679e-07,1.484320847846868e-06
EHF,1.1538425675933854,3.201710496192279e-07,1.5712994766123024e-06
PRRX1,2.0364276706084827,3.2597028765438566e-07,1.5977871392396284e-06
LARP6,1.034522740370137,3.3300533673338904e-07,1.6302595727483684e-06
MMP2,2.159218322190555,3.3874620818292507e-07,1.655814798670288e-06
TRBJ1-6,1.3072002927864221,3.540856625069325e-07,1.7228483903306393e-06
GJA1,1.1362327387143147,3.716603051482994e-07,1.799449241787428e-06
DSG2,1.0203810598192562,4.13924701316396e-07,1.986264470413826e-06
MMP14,1.4472749600134758,4.276892922024024e-07,2.047139111366874e-06
RASSF6,1.3252978623599128,4.7263877483547293e-07,2.2373633010079445e-06
MTCL1,1.230231233785647,5.005361991693618e-07,2.3553947371214258e-06
RASSF8,1.1188799094909336,5.011576219027623e-07,2.357621058824004e-06
P3H2,1.753352709445555,5.123644586531305e-07,2.4060695523634155e-06
STMN2,2.4105000289648872,5.497506909402253e-07,2.569202081950663e-06
PDLIM4,1.0968214486576344,5.51147730127519e-07,2.5737497996789197e-06
GPRIN3,1.1859897603637617,5.547281460492048e-07,2.5889485054604107e-06
FZD7,1.3116469085508022,5.705992410634037e-07,2.653670431815057e-06
AC116351.1,2.2964291350024446,5.793481216869432e-07,2.691209159069447e-06
IL33,1.363512584385954,5.816055899927267e-07,2.7009063565656978e-06
BICC1,1.3705376121192159,6.066035511357741e-07,2.807152852650441e-06
FMO2,1.1378941641789508,6.226169608623746e-07,2.873323173954427e-06
ADAMTS12,1.8409431564463943,6.455571669616761e-07,2.9735661333013045e-06
RASSF2,1.413732745185325,6.709037106493381e-07,3.0840617074846383e-06
CPXM1,1.2525695289420506,6.755397725544545e-07,3.1035780995738755e-06
ITGB4,1.5425384694102338,7.063513130399459e-07,3.2357811158740577e-06
ARNT2,1.2111726837492327,7.758035825893637e-07,3.5315503289921322e-06
DDR1,1.423840584528059,8.279408155267281e-07,3.750626162329116e-06
TC2N,1.3684024632959606,8.828603633111697e-07,3.974496838061272e-06
LAMC2,3.3926824539379594,9.015903607861628e-07,4.049641247543613e-06
FHOD3,1.7561309503547538,9.642920837361305e-07,4.3069364886399355e-06
FAM163B,1.501583907476446,9.996674207106625e-07,4.446198466596724e-06
MITF,1.4147289776522043,1.0361755398365636e-06,4.586757555344477e-06
FOSB,1.8804144936964233,1.0534696271445074e-06,4.652946149788577e-06
FBLN2,2.111567983749889,1.0610759190940796e-06,4.681338553960446e-06
CXCL6,2.392896037367781,1.0634005561637875e-06,4.689484436783571e-06
NRP2,1.123358585287292,1.0845357995818937e-06,4.771597306989699e-06
TNC,2.0825283210177052,1.1263037999518777e-06,4.941685123940118e-06
SGCA,1.4251298661188427,1.147234977930217e-06,5.0251990105204545e-06
SELP,1.1489310069491545,1.1739311456956992e-06,5.120968576283983e-06
SCTR-AS1,1.582418776027948,1.2376594098750634e-06,5.3651033486157645e-06
BHLHE41,1.6767767109753935,1.255704777329262e-06,5.4373966604265725e-06
SGPP2,2.139794568442282,1.3178756642179758e-06,5.675684518669451e-06
BIRC3,1.0898598146394851,1.3437872405842361e-06,5.777885306012176e-06
TNFSF15,1.2587367296790761,1.378305874193202e-06,5.911913161407215e-06
MBOAT2,1.5857397205290704,1.4699882597838265e-06,6.264593398799921e-06
CRHBP,2.375931296330921,1.5052738885940777e-06,6.401239628970428e-06
MLLT3,1.087726740908097,1.605116344659319e-06,6.8021672659836145e-06
SLC7A8,1.3741644116200598,1.7800663468832512e-06,7.481738405620349e-06
SLC2A3,1.0769106634528534,2.2959826659461784e-06,9.450269409083493e-06
HRH1,1.1932214060393842,2.3136315470977886e-06,9.517983196859847e-06
IRAG1,1.5589938825161964,2.430897704551479e-06,9.953455721292864e-06
PKM,1.2158298964663878,2.4752109047693967e-06,1.0104250781256889e-05
ADRA2A,1.908915094802264,2.683592368331005e-06,1.0871147881495835e-05
SLC24A3,2.6567269978697614,2.746888453439575e-06,1.1104918325502412e-05
TRBJ1-5,1.2572788503413943,2.8133580751268797e-06,1.1356307760847354e-05
ANXA3,2.2880928236080673,2.8576175040373486e-06,1.1520336759872607e-05
SLC5A1,3.5925794218999103,3.3040577803600034e-06,1.3136894066183746e-05
AFAP1L2,1.16610925260834,3.4690517222964217e-06,1.3731087195623485e-05
GPC4,1.3614857985670477,3.985342254996753e-06,1.5576838083459156e-05
SAMD11,1.263058105772837,4.562762477236874e-06,1.7668616891146734e-05
ITIH5,3.2371509946161066,4.883044273065067e-06,1.8762598927880663e-05
KLF5,1.1559520690508298,5.062509061715596e-06,1.939591324295761e-05
GPBAR1,1.5059218920951694,5.4326404643079334e-06,2.06992456425309e-05
LIF,1.4670312474824199,5.4980836356572204e-06,2.0918510644424886e-05
EMB,1.0403297148945452,5.914177570828342e-06,2.2314684705883534e-05
LINC02331,1.494225891440307,6.041546612342526e-06,2.27682373043615e-05
KRT80,2.5565726330835634,6.833741679231784e-06,2.5487836164492828e-05
TACSTD2,2.61060238592966,6.939264921579941e-06,2.5857139051438058e-05
MST1R,1.7648150662063553,7.118022002271692e-06,2.6452748642741556e-05
TNFAIP6,3.1223883376630805,7.251684622848551e-06,2.6907744809710684e-05
IGF1R,1.1868448562236595,7.285086873804838e-06,2.7025377116919882e-05
F13A1,1.3842143155601883,7.577863214832394e-06,2.797438775266506e-05
SERPINA3,1.007141806633541,7.6069288086348405e-06,2.806864910447062e-05
MIR6772,1.1057896305934025,7.945318846850778e-06,2.9242599766043063e-05
GLIPR2,1.0005300335298046,8.47021707741775e-06,3.097367625089316e-05
PLA2G4A,1.659248977927296,8.669378293968036e-06,3.162920297473397e-05
SRPX,2.8795328585903857,8.694513000903002e-06,3.171362520935156e-05
CXCL1,2.0764215867962985,9.027894085165963e-06,3.2786704221566865e-05
TFPI2,1.7630815773636197,9.067881440981988e-06,3.2921208268099813e-05
GRAMD1B,1.2163843764590254,9.417093202541798e-06,3.409888691289395e-05
TNFRSF11B,1.3618649252746242,1.0004723711734476e-05,3.605420270272755e-05
CCDC74A,1.2633002766729633,1.0073414436206073e-05,3.6285292517093077e-05
PKHD1,1.7838728747367527,1.1007182597205617e-05,3.9310306574967535e-05
CLIC6,2.457365701275187,1.2694198408582164e-05,4.485151764899358e-05
CHST3,1.4602091943778097,1.4679723530289999e-05,5.133059508083806e-05
FCN3,2.0338539963211195,1.5385340174449755e-05,5.355055121809795e-05
CCDC102B,1.0962173247273577,1.738179555169407e-05,5.9695625977877894e-05
NIBAN1,1.6080943842725768,1.8885731479286654e-05,6.43118225352915e-05
VIPR1,1.1424666479186296,2.0037902784041776e-05,6.782847160458434e-05
SEMA6A,1.0569321976801804,2.0459735406214656e-05,6.907985429422009e-05
CDH6,1.4638277555602173,2.146391337751906e-05,7.214852337683295e-05
ILDR1,1.113276176692574,2.1525764831077992e-05,7.233072263565393e-05
IGFBP5,2.0664724400742513,2.458427692462894e-05,8.164251818425487e-05
CERCAM,1.527092467551459,2.5556362497427082e-05,8.44551156590324e-05
PFKP,1.572178275629855,2.5882387706586773e-05,8.548232598876031e-05
KIT,1.237764528939003,2.790455711708684e-05,9.162793636298542e-05
MRC2,1.8076431194063596,3.278468785418271e-05,0.00010592386140855452
CCN4,1.4979328272500707,3.450594224720853e-05,0.00011110116628754847
OLFM1,1.8113743846425352,3.5035596569314794e-05,0.00011266960207684159
MARCO,2.6107114930976274,3.860151490056946e-05,0.00012304184496880933
NEBL,1.2016510716130173,3.890308353508876e-05,0.00012387888939779064
IL18,1.1504920627536988,4.0041337502942245e-05,0.00012729942018935398
FCN2,2.391190316769395,4.1411602320829694e-05,0.00013128815925247484
DUSP4,1.3934741313133774,4.1971363641625074e-05,0.00013285081570037278
CCDC8,1.2015135186463437,4.341743013097265e-05,0.00013707318160117998
STON1,1.2853753567860846,4.385733396014757e-05,0.0001384070241474103
SLAMF8,1.0421613884387837,4.4154460965520005e-05,0.00013920652747082625
PLAUR,1.4394995248982503,4.477324681112547e-05,0.0001409058664244012
FA2H,2.4906775484901944,4.5183230083197385e-05,0.0001419432066012857
CXCL14,2.2083418169175264,4.78634867856039e-05,0.00014953576767766502
HAVCR2,1.0294850414685282,5.057404696126364e-05,0.00015738018228097557
SPNS3,1.0005271776803617,5.494856192475e-05,0.00016963727721029832
QSOX1,1.0810256874956243,5.65863397620856e-05,0.00017415130889474913
HAND2,2.8514302069637925,5.890830893533215e-05,0.00018073855989886892
TGFB2,1.6920504334176538,6.151507677702843e-05,0.00018793092450603013
PAMR1,1.814801105449632,6.422713026872969e-05,0.00019551023798386195
CXCL5,3.0689769560901277,6.491759778397322e-05,0.0001974985903108207
HUNK,1.4407702876167194,6.748095815702629e-05,0.00020451426327246713
CRP,1.1705231696053364,6.755339048711876e-05,0.00020465574522264912
UNC13D,1.2929597373352653,7.164906234609091e-05,0.0002161574293136195
TMC5,1.7308309692969206,7.326181926597097e-05,0.00022052071180683102
PLCD3,1.079197745857139,7.486835928780299e-05,0.0002247606117542807
CST2,1.1742710952768995,7.508038627633521e-05,0.000225226990045032
DCDC2,1.2296440758992322,7.648422395530577e-05,0.0002290492132617823
SLC52A3,1.1605628850563534,8.235490828657967e-05,0.00024473988074845215
EPHB6,1.281410910515885,8.705639920369014e-05,0.0002571411303589035
TMEM45A,1.0093202317676235,8.765033643023374e-05,0.00025863927007518017
OLFML3,1.0617674976305795,9.051318338256562e-05,0.0002662640680364991
FOSL1,1.6191178245660813,0.0001000518293858089,0.00029192802494802097
UGT1A10,2.458717786226944,0.0001012192038346962,0.0002950092526872627
ADGRG7,1.074547171145483,0.00010423716309519188,0.00030302742219480064
LIPH,1.2500180843159467,0.00010807643446909209,0.0003134433501771005
PLAGL1,1.0098976441538217,0.00010947201012963252,0.0003172592658196242
ADAMDEC1,1.6433983266193128,0.00011130424563110094,0.00032198221811682994
ELOVL7,1.27426176583073,0.00013163643021227737,0.0003745288517369538
GATA3,1.008895476770761,0.00013926366272528943,0.0003943240439403635
C12orf75,1.1478174121226186,0.00014047179387726712,0.0003975324257562824
THSD4,1.0743167875783237,0.00014083604567453819,0.00039835049502534857
SLC35F2,1.1283002688602142,0.00014219495457529067,0.00040176519693011387
CHST4,2.4371049630249666,0.0001470030021313318,0.00041431909962405143
PTX3,1.4892156169839397,0.0001565591386098516,0.0004382994130578023
CX3CR1,1.0489221999647131,0.00016142837304124156,0.000450029010498698
FRZB,1.0297507253623863,0.0001649559162262438,0.0004584961084686783
RAB36,1.706744276525805,0.00016605783732304933,0.0004609141578640112
FZD2,1.3161340687615448,0.00018685078404279254,0.0005121000108869362
IGHA2,1.15305512286534,0.000187836548616939,0.0005143582733531201
LYVE1,1.8627644097987128,0.000190473843094836,0.0005208622414993141
PCDHGB5,1.0414428157728364,0.00019069569672193042,0.0005212895606348764
MEIS3,1.1876948319589993,0.00019251250694145668,0.0005256232927415657
SPINT1,1.4381770746764295,0.0001938836109755097,0.0005290034123011847
COL8A2,2.130387365015311,0.00019978779018486198,0.0005437207487771019
CREB3L1,1.9157673908730932,0.00021529306974981056,0.0005823206800651674
MBOAT4,1.555956221876736,0.00021808069579177405,0.000589158861370843
RERG,1.4578183158202345,0.00022694473232111578,0.0006104083697083682
RAB27B,1.000584114089178,0.0002316010545766669,0.0006215651466403337
PTGFR,1.0331590572026992,0.00025008869788638874,0.0006679116165326
MXRA8,1.832975702578766,0.00028442016651807306,0.0007490296499289744
C2CD4A,1.2523749167778802,0.0002908975599008865,0.0007643153077990895
SH3YL1,1.028127242689023,0.00029583812058204563,0.0007761417337468555
ADRA1A,1.0040051576164681,0.0002962249999647014,0.000777028481755593
SRPX2,1.9643472600200413,0.00030755200227534354,0.0008044835656029061
TIMD4,1.8333079052839876,0.0003094113490374797,0.0008089478296545686
SYT13,2.5167228282449283,0.00031450833240819295,0.0008214631041670449
CDCP1,1.4195247314085393,0.00031469889959073765,0.0008218258103982858
CYS1,2.071620003258083,0.00033152004846073583,0.0008624947120019406
C1orf116,1.4048930708910663,0.00034519837934583123,0.0008944210983017658
DBH,1.0140027483240002,0.00035469287100003895,0.0009169552366636736
SCRN1,1.0592849052881785,0.0003678145351131363,0.0009473071695007152
OLFML2B,1.476332355499123,0.0003692493292416779,0.0009502326918610509
PLXDC1,1.4283921741430055,0.00037505367439515775,0.0009623653281978092
AC010547.2,2.3218384451821548,0.00038422058081919645,0.0009835056928666581"""

print(create_latex_table(data))