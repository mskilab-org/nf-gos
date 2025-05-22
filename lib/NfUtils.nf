
def test_robust_absence(val, test_file = true) {
	def is_absent = ( ! val || (val instanceof Collection && val.empty) )
	def is_present = ! is_absent && val != null
	def is_present_but_fileabsent = is_present && val?.isEmpty()
	if (!test_file) {
		is_present_but_fileabsent = is_absent
	}
	def is_empty = is_absent || is_present_but_fileabsent // Tests if file exists and is nonzero file size
	return is_empty
}

def test_robust_presence(val, test_file = true) {
	def is_absent = ( ! val || (val instanceof Collection && val.empty) )
	def is_present = ! is_absent && val != null
	def is_present_and_filepresent = is_present && ! val?.isEmpty()
	if (!test_file) {
		is_present_and_filepresent = is_present
	}
	return is_present_and_filepresent
}
