package main

import (
	"bytes"
	"fmt"
	"lpp"
	"os"
	//	"sort"
	"flag"
	"strconv"
)

func ParseLine(line []byte) (int, int) {

	line = bytes.TrimSpace(line)
	line_l := bytes.Fields(line)
	depth, err := strconv.Atoi(string(line_l[0]))
	if err != nil {

		panic(fmt.Sprintf(" %s depth not int", (string(line_l[0]))))
	}
	freq, err := strconv.Atoi(string(line_l[1]))
	if err != nil {

		panic(fmt.Sprintf(" %s freq not int", string(line_l[1])))
	}
	return depth, freq
}
func main() {
	kmer := flag.Int("k", 17, "Kmer")
	xls := flag.String("g", "Kmer.xls", "xls File")
	output := flag.String("o", "Kmer.out", "Out File")
	flag.Parse()
	XLS, _ := os.Create(*xls)

	OUT, _ := os.Create(*output)

	all_depth := []int{}
	DATABLOCK := lpp.GetBlockRead("", "\n", false, 1000)
	//	kmer_spieces := make(map[int]float64)
	k_freq := make(map[int]int)
	kmer_individual := make(map[int]int)
	total_kmer_spieces := 0
	total_kmer_number := 0

	error_area := 0
	err_stop := 0
	var maincahe, hetercache, peakmain, valleyheter, peakheter int
	var heter_rate, a_half float64

	valleymain := 0
	start := 0
	//判断错误峰， err_stop以上才用于计算
	for {
		line, err := DATABLOCK.Next()
		if err != nil {
			break
		}
		depth, freq := ParseLine(line)
		kmer_individual[depth] = freq * depth
		if depth > 254 {
			continue
		}
		all_depth = append(all_depth, depth)
		k_freq[depth] = freq

		if depth == 1 && err_stop == 0 {
			error_area += freq * depth
			start = freq
		} else if err_stop == 0 {
			if freq <= start {

				error_area += freq * depth
				start = freq
			} else {
				err_stop = depth - 1
				total_kmer_number += freq * depth
			}

		}
		//如果找到了错误峰,kmer累加，找最大值的peak位置
		if err_stop != 0 {
			total_kmer_spieces += freq
			total_kmer_number += freq * depth
			//			fmt.Println(fmt.Sprintf("%d\t%d", depth, freq*depth))

			if freq*depth > maincahe {
				maincahe = freq * depth
				peakmain = depth
			}
			//计算首个峰的位置
			if peakheter == 0 {
				if depth*freq >= hetercache {
					hetercache = depth * freq
				} else {
					peakheter = depth - 1
				}
			} else if valleyheter == 0 {
				if depth*freq > kmer_individual[depth-1] {
					valleyheter = depth - 1
				}
			}
		}
	}

	//判断最高峰和首个峰是否是一个
	//不存在杂合峰

	if peakheter == peakmain {
		heter_rate = float64(0)
		if valleyheter > 200 {
			valleyheter = 200
		}
		valleymain = valleyheter
	} else {
		//计算主峰的波谷
		for _, dp := range all_depth {
			if dp <= peakmain {
				continue
			}

			if kmer_individual[dp] > kmer_individual[dp-1] {

				valleymain = dp - 1
				break

			}

		}
		if valleymain > 2*peakmain {
			valleymain = 2 * peakmain
		}

		//计算杂合峰的面积
		area_heter := 0
		for i := err_stop; i <= valleyheter; i++ {
			area_heter += kmer_individual[i]

		}

		area_pure := 0
		for i := err_stop; i <= 3*peakmain/2; i++ {
			area_pure += kmer_individual[i]
		}
		a_half = float64(area_heter) / (float64(area_pure) - float64(area_heter))
		//		fmt.Println(a_half)
		heter_rate = a_half / float64(*kmer) / (float64(2) - a_half)
	}
	repeat_kmer := 0
	for i := 2 * peakmain; i <= 255; i++ {
		repeat_kmer += kmer_individual[i]
	}
	repeat_rate := float64(repeat_kmer) / float64(total_kmer_number)
	genome_size := total_kmer_number / (peakmain - 1)
	XLS.WriteString("Depth\tIndiviual Percentage\n")
	for i := err_stop; i <= 3*peakmain; i++ {

		XLS.WriteString(fmt.Sprintf("%d\t%f\n", i, float64(kmer_individual[i])/float64(total_kmer_number)))

	}

	OUT.WriteString("Error Kmer\tTotal Kmer\tUsed Kmer\tPeakmain\tGenomeSize\tRepeatRates\tHeter_Rate\ta[1/2]\n")
	OUT.WriteString(fmt.Sprintf("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n", error_area, total_kmer_number+error_area, total_kmer_number, peakmain, genome_size, repeat_rate, heter_rate, a_half))
}
