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

	"fmt"
	. "lpp"
	"os"

	"github.com/spf13/cobra"
)

func Creeper(taxon_number string) {
	all_end[taxon_number] = ""

	_, ok := taxon_data[taxon_number]
	if ok {
		for k1 := range taxon_data[taxon_number] {
			Creeper(k1)
		}
	}
}

func Taxon_Creep(taxon_number string) {
	all_end = make(map[string]string, 1000000)
	Creeper(taxon_number)
}

var all_end map[string]string = make(map[string]string, 1000000)
var taxon_data map[string]map[string]string = make(map[string]map[string]string, 1000000)

// TaxonCreepCmd represents the TaxonCreep command
var TaxonCreepCmd = &cobra.Command{
	Use:   "TaxonCreep",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	Run: func(cmd *cobra.Command, args []string) {
		file := getFlagString(cmd, "input")
		taxon := getFlagString(cmd, "number")
		output := getFlagString(cmd, "output")
		Taxon_IO := File_Ddict{File_IO: GetBlockRead(*file, "\n", false, 1000000), Header: false}
		taxon_data = Taxon_IO.Read(3, 1)

		Result_File, err := os.Create(*output)
		if err != nil {
			fmt.Printf("%s is not Exist", *output)
			panic("Output Error")
		}
		defer Result_File.Close()
		BufResult := bufio.NewWriter(Result_File)
		defer BufResult.Flush()
		Taxon_Creep(*taxon)
		for k, _ := range all_end {
			BufResult.WriteString(fmt.Sprintf("%s\t%s\n", k, *taxon))
		}
	},
}

func init() {
	rootCmd.AddCommand(TaxonCreepCmd)
	TaxonCreepCmd.PersistentFlags().StringP("input", "i", "", "input Taxon.dump!")
	TaxonCreepCmd.PersistentFlags().StringP("output", "o", "", "Output")
	TaxonCreepCmd.MarkFlagRequired("output")
	TaxonCreepCmd.PersistentFlags().StringP("number", "n", "", "Taxon Number!")

}
