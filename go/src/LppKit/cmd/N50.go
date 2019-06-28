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
	"bytes"

	"fmt"
	. "lpp"
	"os"
	"sort"

	"github.com/spf13/cobra"
)

// N50Cmd represents the N50 command
var N50Cmd = &cobra.Command{
	Use:   "N50",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	Run: func(cmd *cobra.Command, args []string) {
		var all_length int = 0
		//	var length_Ddict map[int]map[string]string = make(map[int]map[string]string)

		file := getFlagString(cmd, "input")
		output := getFlagString(cmd, "output")
		fasta := Fasta{File: *file}
		//		fasta.File = *file

		length_slice := []int{}

		for {
			_, seq, err := fasta.Next()

			//		name := string(data[0])
			seq = bytes.Replace(seq, []byte("\n"), []byte(""), -1)
			seq = bytes.Replace(seq, []byte("\n"), []byte(""), -1)
			seq_length := len(seq)
			fmt.Print(string(seq))

			length_slice = append(length_slice, seq_length)
			all_length = all_length + seq_length
			//		_, ok := length_Ddict[all_length][name]
			//		if !ok {
			//			length_Ddict[all_length] = make(map[string]string)
			//			length_Ddict[all_length][name] = ""
			//		}
			if err != nil {
				break
			}

		}

		sort.Sort(sort.Reverse(sort.IntSlice(length_slice)))
		var N10, N20, N25, N30, N40, N50, N60, N70, N75, N80, N90, Mean int
		var sum_length int = 0

		for _, length := range length_slice {
			sum_length = sum_length + length
			if sum_length >= all_length/10 && N10 == 0 {
				N10 = length

			}
			if sum_length >= all_length/5 && N20 == 0 {
				N20 = length

			}

			if sum_length >= all_length/4 && N25 == 0 {
				N25 = length

			}
			if sum_length >= all_length*3/10 && N30 == 0 {
				N30 = length

			}
			if sum_length >= all_length*2/5 && N40 == 0 {
				N40 = length

			}

			if sum_length >= all_length/2 && N50 == 0 {
				N50 = length

			}
			if sum_length >= all_length*3/5 && N60 == 0 {
				N60 = length

			}
			if sum_length >= all_length*7/10 && N70 == 0 {
				N70 = length

			}
			if sum_length >= all_length*3/4 && N75 == 0 {
				N75 = length

			}
			if sum_length >= all_length*4/5 && N80 == 0 {
				N80 = length

			}
			if sum_length >= all_length*9/10 && N90 == 0 {
				N90 = length

			}

		}

		max := length_slice[0]
		min := length_slice[len(length_slice)-1]
		Mean = sum_length / len(length_slice)
		STATOUT, _ := os.Create(*output + ".stats")
		STATOUT.WriteString("N25\tN50\tN75\tMax\tMin\tMean\tSum.Base\tTotalReads.No\n")
		STATOUT.WriteString(fmt.Sprintf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", N25, N50, N75, max, min, Mean, all_length, len(length_slice)))
		SCOPE, _ := os.Create(*output + ".scope")
		SCOPE.WriteString(fmt.Sprintf("%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n", N10, N20, N30, N40, N50, N60, N70, N80, N90, Mean))
	},
}

func init() {
	rootCmd.AddCommand(N50Cmd)
	N50Cmd.PersistentFlags().StringP("input", "i", "", "Fasta Input")
	N50Cmd.MarkFlagRequired("input")
	N50Cmd.PersistentFlags().StringP("output", "o", "", "Output")
	N50Cmd.MarkFlagRequired("output")
	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// N50Cmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// N50Cmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
