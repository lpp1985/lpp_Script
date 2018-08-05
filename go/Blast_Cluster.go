package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"lpp"
	"os"
	"runtime"
	"sort"
	"strconv"
)

func main() {
	var align_hash map[string]map[string]map[int]string = make(map[string]map[string]map[int]string, 9999999)
	defer func() {
		err := recover()

		if err != nil {

			fmt.Println(err)

		}
	}()
	runtime.GOMAXPROCS(runtime.NumCPU())

	/* Get Option Parser
	 */

	coverage := flag.Float64("c", 0.8, "PAL Result")
	output := flag.String("o", "", "Output")
	fasta := flag.String("f", "", "Fasta Sequence")
	input := flag.String("i", "", "blast_m8")
	flag.Parse()
	Mseq_length := make(map[string]int, 100000000)
	OUTPUT, err := os.Create(*output)

	defer OUTPUT.Close()
	if err != nil {
		panic("Can not Create Result File!!")
	}

	BufOUTPUT := bufio.NewWriterSize(OUTPUT, 9999)
	//	BufOUTPUT.WriteString("TotalBase\tTotalReadsNumber\tQ20%\tQ30%\tN%\tGC%\n")
	//	BufOUTPUT.WriteString(fmt.Sprintf("%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", totalbase, l/4, 100*float64(Q20)/float64(totalbase), 100*float64(Q30)/float64(totalbase), 100*float64(N)/float64(totalbase), 100*float64(GC)/float64(totalbase)))
	defer BufOUTPUT.Flush()
	if err != nil {
		panic("Can not Create Result File!!")
	}

	fasta_io, err := os.Open(*fasta)
	if err != nil {
		panic("Fasta Sequence not exist!!")
	}
	defer fasta_io.Close()
	FASTA_HANDLE := lpp.GetBlockRead(fasta_io, "\n>", false, 10000000)

	for {
		line, err := FASTA_HANDLE.Next()
		if err != nil {
			break
		}
		line = bytes.TrimSuffix(bytes.TrimPrefix(line, []byte(">")), []byte(">"))

		data := bytes.SplitN(line, []byte("\n"), 2)
		name := string(bytes.TrimPrefix(bytes.Fields(data[0])[0], []byte(">")))
		seq := data[1]

		seq = bytes.Replace(seq, []byte("\n"), []byte(""), -1)
		seq_length := len(seq)
		Mseq_length[name] = seq_length

	}
	blast_io, err := os.Open(*input)
	if err != nil {
		panic("blast Sequence not exist!!")
	}
	defer blast_io.Close()
	BLAST_HANDLE := lpp.GetBlockRead(blast_io, "\n", false, 10000000)
	for {
		line, err := BLAST_HANDLE.Next()
		if err != nil {
			break
		}

		line_l := bytes.Split(line, []byte("\t"))
		q_name := string(line_l[0])
		s_name := string(line_l[1])
		//fmt.Println(q_name)
		if q_name == s_name {
			continue
		}

		q_start, _ := strconv.Atoi(string(line_l[6]))
		q_end, _ := strconv.Atoi(string(line_l[7]))
		coor_data := []int{q_start, q_end}
		sort.Ints([]int{q_start, q_end})
		q_start, q_end = coor_data[0], coor_data[1]
		if _, ok1 := align_hash[q_name]; !ok1 {

			align_hash[q_name] = make(map[string]map[int]string, 4000)
			//			align_hash[q_name][s_name] = make(map[int]string, 100000)
		}
		if _, ok2 := align_hash[q_name][s_name]; !ok2 {
			for query, query_hash := range align_hash {
				for sub, sub_hash := range query_hash {
					if float64(len(sub_hash))/float64(Mseq_length[query]) > *coverage {
						fmt.Println("OKAY!!! ", query, sub)
						BufOUTPUT.WriteString(fmt.Sprintf("%s\t%s\t%d\n", query, sub, len(sub_hash)))
					}
				}
			}
			align_hash = make(map[string]map[string]map[int]string, 9999999)
			if float64(q_end-q_start)/float64(Mseq_length[q_name]) < 0.1 {
				continue
			}
			align_hash[q_name] = make(map[string]map[int]string, 4000)
			align_hash[q_name][s_name] = make(map[int]string, 1000000)

		}

		_, ok_start := align_hash[q_name][s_name][q_start]
		_, ok_end := align_hash[q_name][s_name][q_end]
		if ok_start && ok_end {
			continue
		} else {
			for i := q_start; i <= q_end; i++ {
				align_hash[q_name][s_name][i] = ""
			}

		}
	}

	/* Generate Output
	 */

}
