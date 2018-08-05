package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

func getFlagBool(cmd *cobra.Command, flag string) *bool {
	value, err := cmd.Flags().GetBool(flag)
	checkError(err)
	return &value
}

func getFlagString(cmd *cobra.Command, flag string) *string {
	value, err := cmd.Flags().GetString(flag)
	checkError(err)
	return &value
}
func getFlagInt(cmd *cobra.Command, flag string) *int {
	value, err := cmd.Flags().GetInt(flag)
	checkError(err)
	return &value
}
func checkError(err error) {
	if err != nil {
		fmt.Println(err)
		os.Exit(-1)
	}
}
