// Copyright © 2018 NAME HERE <EMAIL ADDRESS>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package cmd

import (
	"bufio"
	"bytes"

	"lpp"
	"os"
	"strconv"

	"github.com/spf13/cobra"
)

// GetSeqCmd represents the GetSeq command
var GetSeqCmd = &cobra.Command{
	Use:   "GetSeq",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	Run: func(cmd *cobra.Command, args []string) {
		var int_hash map[int64]string
		int_hash = make(map[int64]string)
		/* Get Option Parser
		 */
		fasta := getFlagString(cmd, "fasta")

		output := getFlagString(cmd, "output")
		database := getFlagString(cmd, "list")

		listnumber := getFlagInt(cmd, "number")

		seq_number := getFlagBool(cmd, "order")
		header := getFlagBool(cmd, "header")
		start := getFlagBool(cmd, "start")
		exclude := getFlagBool(cmd, "exclude")

		//	defer func() {
		//		if err := recover(); err != nil {
		//			fmt.Println(err)
		//		}
		//	}()
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
	},
}

func init() {
	rootCmd.AddCommand(GetSeqCmd)

	/* Get Option Parser
	 */
	GetSeqCmd.PersistentFlags().StringP("fasta", "f", "", "Fasta Input")

	GetSeqCmd.PersistentFlags().StringP("output", "o", "", "Fasta Output")
	GetSeqCmd.MarkFlagRequired("output")
	GetSeqCmd.PersistentFlags().StringP("list", "d", "", "Need List")
	GetSeqCmd.MarkFlagRequired("list")
	GetSeqCmd.PersistentFlags().IntP("number", "n", 1, "List Number")
	GetSeqCmd.MarkFlagRequired("number")
	GetSeqCmd.PersistentFlags().BoolP("order", "c", false, "Extract Seq from Order")
	GetSeqCmd.PersistentFlags().BoolP("header", "j", false, "Databse has Header")
	GetSeqCmd.PersistentFlags().BoolP("start", "s", false, "Sequence Number start from 1")
	GetSeqCmd.PersistentFlags().BoolP("exclude", "e", false, "Exclude Data from List")

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// GetSeqCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// GetSeqCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
