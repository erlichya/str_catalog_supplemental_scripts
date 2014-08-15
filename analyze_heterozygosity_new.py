import sys
import os
import numpy
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import plotting_functions

def get_counts(lst):
    counts = {}
    for item in lst:
        if item in counts:
            counts[item] = counts[item] + 1.0
        else:
            counts[item] = 1.0
    return counts.values()
    
def calc_heterozyg(freqs):
    val = 1.0
    for freq in freqs:
        val = val - freq*freq
    return val

def read_all_trf_data(src_dir):
    loci  = {}
    files = os.listdir(src_dir)
    for f in files:
        print(f)
        name = f.replace("chr","").replace(".fa", "")
        read_trf_file(loci, name, src_dir+"/"+f)
    return loci
    
def read_trf_file(loci, chrom, input_file):
    data = open(input_file, "r")
    for i in xrange(15):
        data.readline()

    for line in data:
        tokens = line.strip().split()
        patt   = chrom+":"+tokens[0]+"-"+tokens[1]
        period = len(tokens[13])
        score  = int(tokens[7])

        if period < 6:
            loci[patt] = (score, period)
    data.close()
    return loci

def read_reference(input_file):
    data = open(input_file, "r")
    loci = set([])
    for line in data:
        tokens = line.split()
        loci.add(tokens[0])
    data.close()
    return loci

def process_rem_loci(ref_loci, trf_loci, proc_loci, scores, heterozs):
    diffs = [1, 0, -1]
    fail_count    = 0
    success_count = 0
    for locus in ref_loci:
        if locus not in proc_loci:
            should_stop = False

            chrom, other = locus.split(":")
            start, stop  = map(int, other.split("-"))

            for i in diffs:
                if should_stop:
                    break
                for j in diffs:
                    if chrom + ":" + str(start+i) + "-" + str(stop+j) in trf_loci:
                        score, period = trf_loci[chrom + ":" + str(start+i) + "-" + str(stop+j)]
                    
                        if period < 6:
                            heteroz = 0.0
                            heterozs[period-1].append(heteroz)
                            scores[period-1].append(score)

                        success_count = success_count + 1
                        should_stop = True
                        break
            if not should_stop:
                fail_count = fail_count + 1
    print("For remaining loci, TRF entries for "    + str(success_count) + " loci")
    print("For remaining loci, no TRF entries for " + str(fail_count)    + " loci")
    return scores, heterozs


def process_freq_file(trf_loci, input_file):
    data = open(input_file, "r")
    fail_count    = 0
    success_count = 0

    diffs    = [1, 0, -1]
    scores   = [[],[],[],[],[]]
    heterozs = [[],[],[],[],[]]

    proc_loci = set([])

    count_1_total = 0
    count_1 = 0
    vals_1  = []
    count_2_total = 0
    count_2 = 0
    vals_2  = []

    for line in data:
        tokens       = line.replace(";", " ").replace(",", " ").strip().split()
        chrom, ss    = tokens[0].split(":")
        start, stop  = map(int, ss.split("-"))
        allele_lens  = map(lambda x:int(x[0]), map(lambda x:x.split('-'), tokens[1:]))
        freqs        = map(lambda x:float(x[1]), map(lambda x:x.split('-'), tokens[1:]))
        should_stop  = False

        proc_loci.add(chrom+":"+str(start)+"-"+str(stop))

        for i in diffs:
            if should_stop:
                break
            for j in diffs:
                if chrom + ":" + str(start+i) + "-" + str(stop+j) in trf_loci:
                    score, period = trf_loci[chrom + ":" + str(start+i)+"-"+str(stop+j)]
                    
                    if period < 6:
                        heteroz = calc_heterozyg(freqs)
                        heterozs[period-1].append(heteroz)
                        scores[period-1].append(score)

                        if period == 5 and score >= 20 and score < 25:
                            wt_len   = -1
                            max_freq = 0
                            for k in xrange(len(freqs)):
                                if freqs[k] > max_freq:
                                    wt_len   = allele_lens[k]
                                    max_freq = freqs[k]
                            
                            has_multiple = False
                            for k in xrange(len(freqs)):
                                if freqs[k] == max_freq:
                                    continue

                                if (allele_lens[k]-wt_len)%5 == 0:
                                    has_multiple = True
                                    break
                            if has_multiple:
                                count_1 = count_1 + 1
                                vals_1.append(heteroz)
                                print(chrom + "\t" + str(start+i) + "\t" + str(stop+i) + "\t" + str(period) + "\t" + str(score) + "\t" + str(heteroz))

                            count_1_total = count_1_total + 1
                        
                        if period == 5 and score >= 25 and score <= 32:
                            wt_len   = -1
                            max_freq = 0
                            for k in xrange(len(freqs)):
                                if freqs[k] > max_freq:
                                    wt_len   = allele_lens[k]
                                    max_freq = freqs[k]
                            
                            has_multiple = False
                            for k in xrange(len(freqs)):
                                if freqs[k] == max_freq:
                                    continue

                                if (allele_lens[k]-wt_len)%5 == 0:
                                    has_multiple = True
                                    break
                            if has_multiple:
                                count_2 = count_2 + 1
                                vals_2.append(heteroz)
                                print(chrom + "\t" + str(start+i) + "\t" + str(stop+i) + "\t" + str(period) + "\t" + str(score) + "\t" + str(heteroz))

                            count_2_total = count_2_total + 1


                    success_count = success_count + 1
                    should_stop = True
                    break
        if not should_stop:
            fail_count = fail_count + 1
    data.close()

    print("C1 : " + str(count_1) + "\t" + str(count_1_total) + "\t" + str(numpy.sum(map(lambda x: x>0.02, vals_1))))
    print("C2 : " + str(count_2) + "\t" + str(count_2_total) + "\t" + str(numpy.sum(map(lambda x: x>0.02, vals_2))))

   

    #exit(1)

    '''
    pp  = PdfPages("Period_5_plots.pdf")
    fig = plt.figure()
    plotting_functions.plot_distribution(vals_1, fig=fig, title="CDF for P=5 with TRF score from 20-24", xlabel="Heterozygosity", ylabel="CDF", cdf=True)
    pp.savefig(fig)
    fig = plt.figure()
    plotting_functions.plot_distribution(vals_2, fig=fig, title="CDF for P=5 with TRF score from 25-32", xlabel="Heterozygosity", ylabel="CDF", cdf=True)
    pp.savefig(fig)
    pp.close()
    '''

    print("TRF entries for "    + str(success_count) + " loci")
    print("No TRF entries for " + str(fail_count)    + " loci")
    return scores, heterozs, proc_loci


def plot_results(scores, heterozs):
    bins  = numpy.arange(0, 150, 5)
    xvals = [[], [], [], [], []]
    means  = []
    stdevs = []

    cutoffs = [24, 22, 28, 28, 32]
    for i in xrange(5):
        count = 0
        for j in xrange(len(scores[i])):
            if scores[i][j] <= cutoffs[i] and heterozs[i][j] > 0.02:
                count = count + 1
        print(str(i+1)+"\t"+str(count))


    for i in xrange(5):
        c1 = 0
        c2 = 0
        c3 = 0
        c4 = 0
        c5 = 0
        for j in xrange(len(scores[i])):
            c1 = c1 + 1

            if heterozs[i][j] > 0.02:
                if scores[i][j] <= cutoffs[i]:
                    c2 = c2 + 1
                else:
                    c3 = c3 + 1
            else:
                if scores[i][j] <= cutoffs[i]:
                    c4 = c4 + 1
                else:
                    c5 = c5 + 1
        print(i+1, c1, c2, c3, c4, c5)


   

    pp = PdfPages("Heterozygosity_vs_TRF_score_3.pdf")
    for i in xrange(5):
        vals = []
        means.append([])
        stdevs.append([])
        
        for j in xrange(len(bins)):
            vals.append([])
        indices = numpy.digitize(scores[i], bins)
        for j in xrange(len(indices)):
            vals[indices[j]-1].append(heterozs[i][j])
        for j in xrange(len(bins)-1):
            if len(vals[j]) != 0:
                xvals[i].append(0.5*(bins[j]+bins[j+1]))
                means[i].append(numpy.mean(vals[j]))
                stdevs[i].append(numpy.std(vals[j]))
        print(map(len, vals))
        print(sum(map(len, vals)))

    '''
    c1 = 0
    c2 = 0
    c3 = 0
    for j in xrange(len(heterozs[4])):
        if heterozs[i][j] > 0.02:
            if scores[i][j] < 20:
                c1 = c1 + 1
            elif scores[i][j] < 25:
                c2 = c2 + 1
            elif scores[i][j] <= 32:
                c3 = c3 + 1
    print("C's", c1, c2, c3)
    '''

    cutoffs=[24, 22, 28, 28, 32, 34]
    fig = plt.figure()
    for i in xrange(5):
        subplot = fig.add_subplot(230 + i + 1)
        plt.errorbar(xvals[i], means[i], fmt='-o')
        subplot.axvline(x=cutoffs[i], color='r')
        subplot.set_xlabel("TRF Score")
        subplot.set_ylabel("Heterozygosity")
    pp.savefig(fig)
    



    fig = plt.figure()
    for i in xrange(5):
        plt.errorbar(xvals[i], means[i], label=str(i+1), fmt='-o')
    plt.legend()
    plt.xlabel("TRF Score")
    plt.ylabel("Heterozygosity")
    pp.savefig(fig)
    

    for i in xrange(5):
        fig = plt.figure()
        plt.errorbar(xvals[i], means[i], yerr=stdevs[i], label=str(i+1), fmt='-o')
        plt.legend()
        plt.xlabel("TRF Score")
        plt.ylabel("Heterozygosity")
        pp.savefig(fig)
    pp.close()

    
        





def main():
    trf_loci = read_all_trf_data("source_results")
    print("Finished reading TRF information")
    scores, heterozs, proc_loci = process_freq_file(trf_loci, sys.argv[1])
    print("Finished processing frequency file")
    scores, heterozs            = process_rem_loci(read_reference(sys.argv[2]), trf_loci, proc_loci, scores, heterozs)
    print("Finished processing remaining loci")
    plot_results(scores, heterozs)
    


if __name__ == "__main__":
    main()
