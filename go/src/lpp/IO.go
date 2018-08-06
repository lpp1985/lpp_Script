package lpp

import (
	"bufio"
	"bytes"

	"io"

	"os"
)

type Block_Reading struct {
	File     *os.File
	Blocktag string
	Buffer   int
}
type IO struct {
	Io       *bufio.Reader
	BlockTag []byte
	SplitTag byte
}

func GetOuput(name string, buffer int) (*os.File, error) {

	OUTPUTHANDLE, err := os.Create(name)

	OUTPUTIO := bufio.NewWriterSize(OUTPUTHANDLE, buffer)
	defer OUTPUTIO.Flush()

	return OUTPUTHANDLE, err
}

/*
使用说明
FQ2HANDLE, errfq2 := os.Open(*fastq1)
if errfq2 != nil {
	panic(*fastq2 + " Is not Exist!!")
}
FQ2IO := lpp.GetBlockRead(FQ2HANDLE, "\n", false, 10000000)


*/
func (blockreading *Block_Reading) Read() IO {
	BlockIO := IO{}

	raw_file := blockreading.File
	if blockreading.Buffer == 0 {
		blockreading.Buffer = 99999999
	}
	BlockIO.Io = bufio.NewReaderSize(raw_file, blockreading.Buffer)
	if blockreading.Blocktag == "" {
		BlockIO.BlockTag = []byte("\n")
	} else {
		BlockIO.BlockTag = []byte(blockreading.Blocktag)
	}
	BlockIO.SplitTag = byte([]byte(blockreading.Blocktag)[len(blockreading.Blocktag)-1])

	return BlockIO

}

type Fastq struct {
	File string
	IO   IO
}

func (SeqFile *Fastq) Next() (name []byte, seq []byte, name2 []byte, qual []byte, err error) {
	if SeqFile.File == "" {
		SeqFile.IO = GetBlockRead("", "\n>", false, 10000000)
	} else {

		if SeqFile.IO.BlockTag == nil {

			SeqFile.IO = GetBlockRead(SeqFile.File, "\n", false, 10000000)

		}
	}

	name, err = SeqFile.IO.Next()
	seq, err = SeqFile.IO.Next()
	name2, err = SeqFile.IO.Next()
	qual, err = SeqFile.IO.Next()
	return name, seq, name2, qual, err

}

type Fasta struct {
	File string
	IO   IO
}

//func (fa *Fasta) Read() *Fasta {
//	FastaIO := GetBlockRead(fa.File, "\n>", false, 10000000)
//	fa.IO = FastaIO

//	return fa
//}

func (SeqFile *Fasta) Next() (name []byte, seq []byte, err error) {
	if SeqFile.File == "" {
		SeqFile.IO = GetBlockRead("", "\n>", false, 10000000)
	} else {

		if SeqFile.IO.BlockTag == nil {

			SeqFile.IO = GetBlockRead(SeqFile.File, "\n>", false, 10000000)

		}
	}
	Data_block, err := SeqFile.IO.Next()

	Data_block = bytes.TrimSuffix(Data_block, []byte(">"))
	if string(Data_block[0]) != ">" {

		Data_block = append([]byte(">"), Data_block...)
	}

	Data_list := bytes.SplitN(Data_block, []byte("\n"), 2)
	name = append(Data_list[0], '\n')
	seq = Data_list[1]

	return name, seq, err
}
func GetBlockRead(filehandle string, blocktag string, header bool, buffer int) IO {
	BR := new(Block_Reading)
	BR.Blocktag = blocktag
	BR.Buffer = buffer
	FILE := os.Stdin
	if filehandle != "" {

		c_FILE, errfile := os.Open(filehandle)
		if errfile != nil {
			panic(filehandle + " Is not Exist!!")
		}
		FILE = c_FILE
	}
	BR.File = FILE
	Result_IO := BR.Read()

	if header {
		Result_IO.Next()
	}
	return Result_IO
}
func (Reader IO) Next() ([]byte, error) {

	var out_tag []byte
	var status error

	for {
		line, err := Reader.Io.ReadSlice(Reader.SplitTag)
		status = err
		out_tag = append(out_tag, line...)
		if err == nil {

			if len(Reader.BlockTag) > 1 {

				if len(out_tag) >= len(Reader.BlockTag) && bytes.Equal(out_tag[(len(out_tag)-len(Reader.BlockTag)):], Reader.BlockTag) {

					break
				}

			} else {
				break
			}

		} else if err == io.EOF {
			break
		} else if err != bufio.ErrBufferFull {
			break
		}

	}

	return out_tag, status

}
