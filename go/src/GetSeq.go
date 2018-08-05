// GetSeq
package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"lpp"
	"os"
	"strconv"
)

func main() {

	var int_hash map[int64]string
	int_hash = make(map[int64]string)
	/* Get Option Parser
	 */
	fasta := flag.String("f", "", "Fasta Input")
	output := flag.String("o", "", "Fasta Output")
	database := flag.String("d", "", "Need List")
	listnumber := flag.Int("n", 1, "List Number")
	seq_number := flag.Bool("c", false, "Extract Seq from Order")
	header := flag.Bool("j", false, "Databse has Header")
	start := flag.Bool("s", false, "Sequence Number start from 1")
	exclude := flag.Bool("e", false, "Exclude Data from List")
	flag.Parse()
	defer func() {
		if err := recover(); err != nil {
			fmt.Println(err)
		}
	}()
	/* Generate Output
	 */
	Result_File, err := os.Create(*output)
	if err != nil {
		panic("Can not Create Output File!!")
	}
	defer Result_File.Close()
	BufResult := bufio.NewWriter(Result_File)
	defer BufResult.Flush()
	/*
		open Fasta Input

	*/
	seqIO := lpp.Fasta{File: *fasta}

	/*
		Prepeare Data Input Database
	*/

	dataIO := lpp.GetBlockRead(*database, "\n", *header, 1000000)
	if err != nil {
		panic("Database error")
	}

	/* need_hash Generate
	 */
	raw_hash := new(lpp.File_dict)
	raw_hash.File_IO = dataIO
	raw_hash.Header = false
	need_hash := raw_hash.Read(*listnumber, *listnumber)

	if *seq_number {
		for t, _ := range need_hash {
			number, err := strconv.Atoi(t)
			if err == nil {
				number2 := int64(number)
				int_hash[number2] = ""
			} else {
				panic("Input Database not Integer")
			}

		}
	}

	var i int64 = 0
	if *start {
		i = 1
	}

	for {
		//fmt.Println(i)
		title, sequence, err := seqIO.Next()

		var ok bool = false
		if *seq_number {

			_, has := int_hash[i]
			ok = has
			if *exclude {
				ok = !has
			}

		} else {

			name := bytes.Fields(title)[0][1:]

			_, has := need_hash[string(name)]

			ok = has
			if *exclude {
				ok = !has
			}

		}
		if ok {

			BufResult.Write(title)
			BufResult.Write(sequence)

		}
		if err != nil {
			break
		}
		i++

	}

}
