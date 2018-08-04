#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/7/24
"""

from lpp import * 



if __name__ == '__main__':
	all_file = sys.argv[1:-1]
	data_group = {}
	for e_f in all_file:
		sample = e_f.split(".")[0]
		if sample not in data_group:
			data_group[sample] = [e_f]
		else:
			data_group[sample].append(e_f)
			
	print(data_group)
	END = open(sys.argv[-1], 'w')
	END.write("Sample\tRawBase\tRawReadsNumber\tRaw_Q20%\tRaw_Q30%\tRAW_N%\tRaw_GC%\tFilteredBase\tFilteredReadsNumber\tFiltered_Q20%\tFiltered_Q30%\tRAW_N%\tFiltered_GC%\tQC_DataPerc%\tQC_ReadsPerc%\n")
	for sample in data_group:
		Filterdata, Rawdata = sorted(data_group[sample])
		RAW =  open(Rawdata, 'rU')
		RAW.next()
		line = RAW.next()
		line_l = line.split("\t")
		raw_base = int(line_l[0])
		raw_number =   int(line_l[1])
		END.write(sample + '\t' +line.strip() + '\t')
		RAW =  open(Filterdata, 'rU')
		RAW.next()
		line = RAW.next()
		line_l = line.split("\t")
		qc_base = int(line_l[0])
		qc_number =   int(line_l[1])
		base_perc = "%.2f" % (100.0 * qc_base / raw_base)
		number_perc = "%.2f" % (100.0 * qc_number / raw_number)
		END.write( line.strip() + '\t' + base_perc + '\t' + number_perc + '\n')
