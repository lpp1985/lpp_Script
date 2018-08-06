package lpp

func RevComplement(char []byte) []byte {
	var complement = [256]uint8{
		'A': 'T', 'a': 'T',
		'C': 'G', 'c': 'G',
		'G': 'C', 'g': 'C',
		'T': 'A', 't': 'A',
		'U': 'A', 'u': 'A',
		'M': 'K', 'm': 'K',
		'R': 'Y', 'r': 'Y',
		'W': 'W', 'w': 'W',
		'S': 'S', 's': 'S',
		'Y': 'R', 'y': 'R',
		'K': 'M', 'k': 'M',
		'V': 'B', 'v': 'B',
		'H': 'D', 'h': 'D',
		'D': 'H', 'd': 'H',
		'B': 'V', 'b': 'V',
		'N': 'N', 'n': 'N',
	}
	L := len(char)
	new_base := make([]byte, L)

	for _, base := range char {
		L--
		new_base[L] = complement[base]

	}
	return new_base
}
