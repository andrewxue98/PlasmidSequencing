# from https://github.com/ponnhide/PySanger

from Bio import SeqIO
from Bio import pairwise2
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 

from codon_table import dna_codons
from colors import rasmol

matplotlib.rcParams['font.family']       = 'sans-serif'
matplotlib.rcParams['font.sans-serif']   = ["Arial","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['figure.figsize']    = [3,3]
matplotlib.rcParams['font.size']         = 10
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 1.0 
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0

_atgc_dict = {0:"A", 1:"T", 2:"G", 3:"C"}

def abi_to_dict(filename):
    record   = SeqIO.read(filename,'abi')
    abi_data = {"conf":[],
                "channel":{"A":[],
                           "T":[],
                           "G":[],
                           "C":[],
                          },
                "_channel":{"A":[],
                            "T":[],
                            "G":[],
                            "C":[],
                          }
                }
    for i, (pos, conf) in enumerate(zip(record.annotations['abif_raw']["PLOC1"], record.annotations['abif_raw']["PCON1"])):
        if pos > 4 and pos < len(record.annotations['abif_raw']["DATA9"])-5: 
            abi_data["conf"].append(conf)
            abi_data["channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])

            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-5])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+5])
            
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-5])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+5])
            
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-5])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+5])

            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-5])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+5])

    return abi_data 

def generate_consensusseq(abidata, conf_thresh = 10):
    consensus_seq = "" 
    
    for values in zip(abidata["channel"]["A"], abidata["channel"]["T"], abidata["channel"]["G"], abidata["channel"]["C"]):
        consensus_seq += _atgc_dict[values.index(max(values))]

    conf_consensus_seq = ""
    for conf, base in zip(abidata["conf"], consensus_seq):
        if conf > conf_thresh:
            conf_consensus_seq += base
        else:
            conf_consensus_seq += "N"
     
    return (conf_consensus_seq, conf_consensus_seq.translate(str.maketrans("ATGCN","TACGN"))[::-1]) 

def generate_pwm(abidata):
    pwm = {"A":[], "T":[], "G":[], "C":[]} 
    for values in zip(abidata["channel"]["A"], abidata["channel"]["T"], abidata["channel"]["G"], abidata["channel"]["C"]):
        v = 100000 / (sum(values)+1) 
        new_values = (v*values[0], v*values[1], v*values[2], v*values[3])
        new_values = list(map(int, new_values))
        
        while sum(new_values) < 100000:
            for i in range(len(new_values)):
                new_values[i] += 1
                if sum(new_values) == 100000:
                    break 
        
        pwm["A"].append(new_values[0])
        pwm["T"].append(new_values[1])
        pwm["G"].append(new_values[2])
        pwm["C"].append(new_values[3])
    
    pwm=pd.DataFrame(pwm)
    return pwm 

def _colorbar(ax, ref, matches=None, char=True, fontsize=8, width = 0.96):
    bars = ax.bar([item + 0.5 for item in range(len(ref))], [0.9] * (len(ref)), width=width, edgecolor="#BBBBBB", linewidth=0.3, align="center",bottom=0.05)
    ax.set_xlim(0,len(ref))
    ax.set_ylim(0,1.00)
    p = 0
    if matches is None:
        for bar, c in zip(bars,ref):
            bar.set_facecolor("w")
            if char == True:
                ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
            p += 1
    
    else:
        for m, bar, c in zip(matches, bars,ref):
            bar.set_edgecolor("#BBBBBB")
            bar.set_linewidth(0.5)
            if m == 1:
                bar.set_facecolor("#FFFFFF")
            if m == 0:
                bar.set_facecolor("#FF0000")
            if m == -1:
                bar.set_facecolor("#C00000")
            if m == -2:
                bar.set_facecolor("#FFC000")
                    
            if char == True:
                ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
            p += 1

    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #ax.patch.set_alpha(0.0)
    return ax

def _aabar(ax, ref, mismatches=None, char=True, fontsize=8, width = 0.96):
    last_bar_pos = 0
    bar_pos = []
    widths = []
    total = 0

    for item in ref:
        if item == '*':
            last_bar_pos += 0.5/3
            bar_pos.append(last_bar_pos)
            last_bar_pos += 0.5/3
            widths.append(width/3)
        else:
            last_bar_pos += 0.5
            bar_pos.append(last_bar_pos)
            last_bar_pos += 0.5
            widths.append(width)

    total = last_bar_pos
    
    bars = ax.bar(bar_pos, [0.9] * len(ref), width=widths, edgecolor="#BBBBBB", linewidth=0.3, align="center",bottom=0.05)
    #bars = ax.bar([item + 0.5 for item in list(range(len(ref)))], [0.9] * len(ref), width=width, edgecolor="#BBBBBB", linewidth=0.3, align="center",bottom=0.05)
    ax.set_xlim(0, total)
    ax.set_ylim(0, 1.00)

    p = 0
    if mismatches is None:
        for bar, c in zip(bars, ref):
            if c == '*':
                bar.set_facecolor("w")
                p += 1/3
            else:
                if c in rasmol:
                    bar.set_facecolor(rasmol[c])
                else:
                    bar.set_facecolor(rasmol["Other"])

                if char == True:
                    ax.text(p+ 0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
                    
                p += 1

    else: 
        for bar, c, m in zip(bars, ref, mismatches):
            if c == '*':
                bar.set_facecolor("w")
                p+=1/3
            else:
                if m == -1:
                    bar.set_facecolor("#f5f5f5")
                elif m == 0:
                    bar.set_facecolor("#f5f5f5")
                elif m == 1:
                    bar.set_facecolor("#408040")
                elif m == 2:
                    bar.set_facecolor("#808040")
                elif m == 3:
                    bar.set_facecolor("#804040")
                elif m == 4:
                    bar.set_facecolor("#802020")
                elif m == 5:
                    bar.set_facecolor("#801010")
                else:
                    bar.set_facecolor("w")
                if char == True:
                    ax.text(p+ 0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
                p += 1
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    return ax

def visualize(title, abidata, template=None, strand=1, fig=None, region="read", translation_limits = None, template_limits = None):  
    '''template limits should be inclusive'''
    avalues = abidata["_channel"]["A"]
    tvalues = abidata["_channel"]["T"]
    gvalues = abidata["_channel"]["G"]
    cvalues = abidata["_channel"]["C"]
    consensus_seq_set = generate_consensusseq(abidata)

    axes_width = 1 / 60

    if strand == 1: 
        subject = consensus_seq_set[0]

    if strand == -1:
        subject = consensus_seq_set[1] 
        avalues, tvalues, gvalues, cvalues = tvalues[::-1], avalues[::-1], cvalues[::-1], gvalues[::-1]

    if template_limits is not None:
        template_start = max(0, template_limits[0]-1) #adjust for 0-based indexing
        if template_limits[1] is None:
            template_end = len(template)
        else:
            template_end = min(len(template), template_limits[1])

        if template_start >= template_end:
            raise ValueError("Invalid template limits specified! Template start point should be less than template end point.")
    else:
        template_start = 0
        template_end = len(template)

    if translation_limits is not None:
        translation_start = max(0, translation_limits[0]-1)
        if translation_limits[1] is None:
            translation_end = len(template)
        else:
            translation_end = min(len(template), translation_limits[1])

    template = template[template_start:template_end]

    alignments = pairwise2.align.globalms(template, subject, 2, 0, -10, -1, penalize_end_gaps=False)
    atemplate = alignments[0][0]
    asubject = alignments[0][1]

    new_avalues = []
    new_tvalues = [] 
    new_gvalues = []
    new_cvalues = []
    new_conf = []

    #trim off leading and trailing gaps, align front and end
    ts = len(atemplate) - len(atemplate.lstrip('-'))
    te = len(atemplate) - len(atemplate.rstrip('-'))

    ss = len(asubject) - len(asubject.lstrip('-'))
    se = len(asubject) - len(asubject.rstrip('-'))

    #modify to proper indices
    te = len(atemplate) - te 
    se = len(asubject)  - se

    #modify actg values to include gaps
    pos = 0
    for s in asubject:
        if s == "-":
            new_avalues.extend([0,0,0,0,0]) 
            new_tvalues.extend([0,0,0,0,0]) 
            new_gvalues.extend([0,0,0,0,0]) 
            new_cvalues.extend([0,0,0,0,0]) 
            new_conf.append(0)
        else:
            new_avalues.extend(avalues[5*pos:5*pos+5])
            new_tvalues.extend(tvalues[5*pos:5*pos+5])
            new_gvalues.extend(gvalues[5*pos:5*pos+5])
            new_cvalues.extend(cvalues[5*pos:5*pos+5])
            new_conf.append(abidata["conf"][pos])
            pos += 1

    #find mismatches #1 is match, 0 is mismatch, -1 is gap, -2 is N
    matches = []
    for t,s in zip(atemplate.upper(), asubject.upper()): 
        if t==s:
            matches.append(1)
        elif t == '-' or s == '-':
            matches.append(-1)
        elif t.upper() == 'N' or s.upper() == 'N':
            matches.append(-2)
        else:
            matches.append(0)

    if region == "template":
        start, end = ts, te
    elif region == 'read':
        start, end = ss, se
    else:
        raise ValueError("Invalid region specified!")

    asubject  = asubject[start:end]
    atemplate = atemplate[start:end]
    matches   = matches[start:end] 

    if len(asubject) != len(atemplate):
        raise ValueError("Subject and template lengths do not match after truncation!")
    
    length = len(asubject)

    if region == 'read':
        t_start_base = max(template_start, template_start+ss-ts)
        s_start_base = 0
    elif region == 'template':
        t_start_base = template_start
        s_start_base = max(0, ts-ss)

    remove_start_num = len(atemplate) - len(atemplate.lstrip('-'))
    at = atemplate.lstrip('-')
    asub = asubject[remove_start_num:]

    starting_codon = 1
    starting_base = t_start_base + 1 #return to 1 index

    def remove_index(inp_str, remove):
        char_removed = 0
        for idx in range(len(inp_str)):
            if inp_str[idx] != '-':
                char_removed += 1
            
            if char_removed == remove:
                return idx + 1


    if translation_start > t_start_base:
        remove = translation_start - t_start_base #need to cut to remove this many bases from the start of the template
    else:
        remove = (translation_start - t_start_base) % 3 #need to cut to remove this many bases from the start of the template
        starting_codon = starting_codon - (remove // 3) #rounds up

    starting_base += remove

    at = at[remove_index(at, remove):]
    asub = asub[remove_index(at, remove):]

    temp_chars = len(at.rstrip('-')) - at.rstrip('-').count('-') #number of characters in template minus number of gaps
    temp_chars = temp_chars // 3 * 3 #round down to nearest multiple of 3

    at = at[:remove_index(at, temp_chars)]
    asub = asub[:remove_index(at, temp_chars)]

    #cluster atemplate into 3 letter chunks from template
    at_codons = []
    asub_codons = []
    template_nt_pos = []
    codon_num = []

    current_codon = starting_codon
    current_nt_pos = starting_base
    curr_template_codon = ''
    curr_subject_codon = ''
    n_template_nt = 0
    curr_codon_pos = []

    for tc, sc in zip(at, asub):
        curr_template_codon += tc
        curr_subject_codon += sc
        curr_codon_pos.append(current_nt_pos)
        
        if tc != '-':
            n_template_nt += 1
            current_nt_pos += 1
        
        if n_template_nt == 3:
            at_codons.append(curr_template_codon)
            asub_codons.append(curr_subject_codon)
            template_nt_pos.append(curr_codon_pos)
            codon_num.append(current_codon)
            
            curr_template_codon = ''
            curr_subject_codon = ''
            curr_codon_pos = []
            n_template_nt = 0
            current_codon += 1    

    if n_template_nt != 0:
        print(n_template_nt)
        raise ValueError('Something has gone wrong :(')

    insertions = []
    deletions = []
    nt_mismatches = []
    silent_mutations = []
    missense_mutations = []
    frameshift_mutations = []

    for at_codon, as_codon, codon_pos, curr_codon_num in zip(at_codons, asub_codons, template_nt_pos, codon_num):
        for at_codon_nt, as_codon_nt, nt_pos in zip(at_codon, as_codon, codon_pos):
            if at_codon_nt == '-':
                insertions.append(f"{nt_pos}(ins{as_codon_nt.upper()})")
            elif as_codon_nt == '-':
                deletions.append(f"{nt_pos}(del{at_codon_nt.upper()})")
            elif at_codon_nt != as_codon_nt:
                nt_mismatches.append(f"{nt_pos}({at_codon_nt.upper()}->{as_codon_nt.upper()})")

        at_codon_s = at_codon.strip('-').upper()
        as_codon_s = as_codon.strip('-').upper()

        if len(as_codon_s) != 3 or len(at_codon_s) != 3:
            frameshift_mutations.append(f"{curr_codon_num}")
        else:
            #both codons are 3 nt long
            at_aa = translate([at_codon_s])[0]
            as_aa = translate([as_codon_s])[0]

            if at_codon_s != as_codon_s:
                if at_aa == as_aa:
                    silent_mutations.append(f"{curr_codon_num}")
                else:
                    missense_mutations.append(f"{at_aa}{curr_codon_num}{as_aa}")

    insertions = "; ".join(insertions)
    deletions = "; ".join(deletions)
    nt_mismatches = "; ".join(nt_mismatches)
    silent_mutations = "; ".join(silent_mutations)
    missense_mutations = "; ".join(missense_mutations)
    frameshift_mutations = "; ".join(frameshift_mutations)

    stats = pd.DataFrame({'id': title, 'alignment_score': alignments[0].score, 'nt_insertions' : [insertions], 'nt_deletions': [deletions], 'nt_mismatches': [nt_mismatches], 'silent_mutations': [silent_mutations], 'missense_mutations': [missense_mutations], 'gap_mutations': [frameshift_mutations]})

    translation_start = max((translation_start - t_start_base) % 3, translation_start + ts - t_start_base) #make sure translation start stays in right reading frame
    translation_end = min(len(atemplate), translation_end + ts - t_start_base)

    avalues = new_avalues[5*start:5*end]
    tvalues = new_tvalues[5*start:5*end]
    gvalues = new_gvalues[5*start:5*end]
    cvalues = new_cvalues[5*start:5*end]
    conf = new_conf[start:end]

    if fig is None:
        fig = plt.figure(figsize=(200, 3.5))

    #dynamically scale figsize according to 

    atemplate_codons = make_codons(atemplate, translation_start=translation_start, translation_end=translation_end)
    atemplate_aa = translate(atemplate_codons)
    asubject_codons = make_codons(asubject, translation_start=translation_start, translation_end=translation_end)
    asubject_aa = translate(asubject_codons)

    aa_mismatches = []
    for t_aa, s_aa, t_c, s_c in zip(atemplate_aa, asubject_aa, atemplate_codons, asubject_codons):
        if t_aa == s_aa:
            if t_aa == '-' or t_aa == 'X':
                aa_mismatches.append(-1) #gap or unknown, not real match or mismatch
            elif t_c == s_c:
                aa_mismatches.append(0) #no mutation
            else:
                aa_mismatches.append(1) #silent mutation
        elif t_aa == 'X' or s_aa == 'X':
            aa_mismatches.append(2)
        elif t_aa == '-' or s_aa == '-':
            aa_mismatches.append(4)
        elif t_aa == 'FS' or s_aa == 'FS':
            aa_mismatches.append(5)
        else:
            aa_mismatches.append(3)

    ax  = fig.add_axes([0, 0.5, axes_width * length/20,  0.3])
    axs = fig.add_axes([0, 0.34, axes_width * length/20, 0.08])
    axst = fig.add_axes([0, 0.42, axes_width * length/20, 0.08])
    axst = _aabar(axst, asubject_aa, aa_mismatches, char=True)
    axst.set_ylabel("AA Matches", rotation=0, va="center", ha="right", fontsize = 'small') 

    ax.bar(list(range(length)), conf, width=1.0, edgecolor="#BBBBBB", linewidth=0.5, facecolor="#F7F7FF", align="edge", zorder=0)
    ax.set_xlim(0,length)

    ax2 = ax.twinx()
    positions = list(map(lambda x: (x+0.5)/5, list(range(len(tvalues)))))
    ax2.plot(positions, tvalues, color="#FC58FE", lw=1, zorder=1)  
    ax2.plot(positions, avalues, color="#33CC33", lw=1, zorder=1) 
    ax2.plot(positions, gvalues, color="#303030", lw=1, zorder=1)
    ax2.plot(positions, cvalues, color="#395CC5", lw=1, zorder=1) 

    all_values = list(tvalues) + list(avalues) + list(gvalues) + list(cvalues)
    ax2.set_ylim(min(all_values), 1.01*max(all_values))

    ax.set_ylabel("Quality")
    ax2.set_ylabel("Peak Height") 

    ax.spines["top"].set_visible(False) 
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False) 
    ax2.spines["top"].set_visible(False) 
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["left"].set_visible(False)

    axs = _colorbar(axs, asubject, matches=matches, char=True)
    tick_space = min(max(length // 20, 1), 10)

    

    ticks, ticklabels = [], []

    base_num = 0
    pos_num = 0
    for s in asubject:
        if s != '-':
            if(base_num % tick_space == 0):
                ticks.append(pos_num + 0.5)
                ticklabels.append(s_start_base + base_num + 1)
            base_num += 1
        pos_num += 1
        
    axs.set_xticks(ticks) 
    axs.set_xticklabels(ticklabels) 

    ax.set_xticks([]) 
    ax.set_xticklabels([]) 
    axs.set_ylabel("Read", rotation=0, va="center", ha="right", fontsize = 'small') 

    axt = fig.add_axes([0, 0.88, axes_width * length/20, 0.08])
    #axat = fig.add_axes([0, 0.88, axes_width * length/20, 0.08])
    axtt = fig.add_axes([0, 0.8, axes_width * length/20, 0.08])
    axtt = _aabar(axtt, atemplate_aa, char=True)
    axtt.set_ylabel("Template AA", rotation=0, va="center", ha="right", fontsize = 'small') 

    axt = _colorbar(axt, atemplate, matches=matches, char=True)
    ticks, ticklabels = [],[]

    base_num = 0
    pos_num = 0
    for t in atemplate:
        if t != '-':
            if(base_num % tick_space == 0):
                ticks.append(pos_num + 0.5)
                ticklabels.append(t_start_base + base_num + 1)
            base_num += 1
        pos_num += 1

    axt.xaxis.tick_top() 
    axt.set_xticks(ticks) 
    axt.set_xticklabels(ticklabels) 

    axt.set_ylabel("Template", rotation=0, va="center", ha="right", fontsize = 'small') 
    fig.suptitle(title, x = 0, y = 1, fontsize=10)

    return fig, stats

def translate(codons):
    protein = []
    for codon in codons:
        if codon == '*':
            protein.append('*') #temporary filler character
        else:
            if codon in dna_codons:
                protein.append(dna_codons[codon])
            elif codon == "---" or len(codon) < 3:
                protein.append("-")
            elif codon.count('-') > 0:
                protein.append("FS")
            elif codon.count('N') > 0:
                protein.append('X')
            else:
                raise ValueError(f'Unrecognized codon: {codon}')
    return protein

def make_codons(seq, translation_start = 0, translation_end = None, fill_char = '*'):
    seq = seq.upper()
    codons = []
    if translation_end == None:
        translation_end = len(seq)
    translation_end = translation_end - (translation_end-translation_start) % 3

    for h in range(translation_start):
        codons.append(fill_char)
    for i in range(translation_start, translation_end, 3):
        codon = seq[i:i+3]
        codons.append(codon)
    for j in range(translation_end, len(seq)):
        codons.append(fill_char)
    return codons

