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

	"fmt"
	"lpp"
	"os"

	"github.com/spf13/cobra"
)

func QC(datahandle *string, qc_number byte, c chan [6]int) {
	l := 0
	Q20 := 0
	Q30 := 0
	GC := 0
	N := 0
	totalbase := 0

	dataIO := lpp.GetBlockRead(*datahandle, "\n", false, 10000000)
	for {

		line, err := dataIO.Next()
		if err != nil {
			break
		}
		l += 1
		if l%4 == 0 {
			line_cont := bytes.Split(line, []byte("\n"))[0]

			for i := 0; i < len(line_cont); i++ {
				if line_cont[i]-qc_number > 20 {
					Q20 += 1
					if line_cont[i]-qc_number > 30 {
						Q30 += 1
					}
				}

			}
		} else if l%2 == 0 {

			line_cont := bytes.Split(line, []byte("\n"))[0]
			totalbase += len(line_cont)
			for i := 0; i < len(line_cont); i++ {
				if string(line_cont[i]) == "N" {
					N += 1

				} else if string(line_cont[i]) == "G" || string(line_cont[i]) == "C" {

					GC += 1
				}
			}
		}

	}
	c <- [6]int{l, totalbase, Q20, Q30, GC, N}
}

// QualityCatCmd represents the QualityCat command
var QualityCatCmd = &cobra.Command{
	Use:   "QualityCat",
	Short: "Staistics quality score for a fastq file!!",
	Long:  `Staistics quality score for a fastq file!!`,
	Run: func(cmd *cobra.Command, args []string) {
		result_channel := make(chan [6]int)
		read1 := getFlagString(cmd, "1")
		read2 := getFlagString(cmd, "2")
		output := getFlagString(cmd, "output")
		quality := getFlagInt(cmd, "number")
		qc_number := byte(*quality)

		l := 0
		Q20 := 0
		Q30 := 0
		GC := 0
		N := 0
		totalbase := 0
		go QC(read1, qc_number, result_channel)
		go QC(read2, qc_number, result_channel)
		data1, data2 := <-result_channel, <-result_channel
		l = data1[0] + data2[0]
		totalbase = data1[1] + data2[1]
		Q20 = data1[2] + data2[2]
		Q30 = data1[3] + data2[3]
		GC = data1[4] + data2[4]
		N = data1[5] + data2[5]
		//	fmt.Println(data1, data2)
		/* Generate Output
		 */
		OUTPUT, err := os.Create(*output + ".stats")

		if err != nil {
			panic("Can not Create Result File!!")
		}

		BufOUTPUT := bufio.NewWriterSize(OUTPUT, 9999)
		BufOUTPUT.WriteString("TotalBase\tTotalReadsNumber\tQ20%\tQ30%\tN%\tGC%\n")
		BufOUTPUT.WriteString(fmt.Sprintf("%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", totalbase, l/4, 100*float64(Q20)/float64(totalbase), 100*float64(Q30)/float64(totalbase), 100*float64(N)/float64(totalbase), 100*float64(GC)/float64(totalbase)))
		BufOUTPUT.Flush()
	},
}

func init() {
	rootCmd.AddCommand(QualityCatCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// QualityCatCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// QualityCatCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")

	/* Get Option Parser
	 */
	QualityCatCmd.PersistentFlags().StringP("1", "1", "", "Read1")
	QualityCatCmd.PersistentFlags().StringP("2", "2", "", "Read2")
	QualityCatCmd.PersistentFlags().StringP("output", "o", "", "Output")
	QualityCatCmd.MarkFlagRequired("output")
	QualityCatCmd.PersistentFlags().IntP("number", "n", 33, "Qualit Score")

}
