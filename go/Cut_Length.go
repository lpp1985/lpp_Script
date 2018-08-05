package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"lpp"
	"os"
	"runtime"
)

func main() {
	defer func() {
		err := recover()

		if err != nil {

			fmt.Println(err)

		}
	}()
	runtime.GOMAXPROCS(runtime.NumCPU())

	/* Get Option Parser
	 */

	length := flag.Int("l", 0, "PAL Result")
	output := flag.String("o", "", "Output")
	other := flag.String("r", "", "Other")
	fasta := flag.String("f", "", "Fasta Sequence")
	flag.Parse()
	OTHER, err := os.Create(*other)
	defer OTHER.Close()
	BufOTHER := bufio.NewWriterSize(OTHER, 9999)
	defer BufOTHER.Flush()

	if err != nil {
		panic("Can not Create Other File!!")
	}
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
		line = bytes.TrimSuffix(bytes.TrimPrefix(line, []byte(">")), []byte(">"))

		data := bytes.SplitN(line, []byte("\n"), 2)
		name := string(bytes.TrimPrefix(bytes.Fields(data[0])[0], []byte(">")))
		seq := data[1]

		seq = bytes.Replace(seq, []byte("\n"), []byte(""), -1)
		seq_length := len(seq)
		result := []byte(">")
		result = append(result, name...)
		result = append(result, []byte("\n")...)
		result = append(result, seq...)
		result = append(result, []byte("\n")...)
		if seq_length > *length {

			BufOUTPUT.Write(result)
		} else {
			BufOTHER.Write(result)
		}

		if err != nil {
			break
		}
	}

	//	fmt.Println(Mseq_length)
	/* Generate Output
	 */

}
