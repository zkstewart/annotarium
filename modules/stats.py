#! python3

def N50(numlist): 
  """ 
  Abstract: Returns the N50 value of the passed list of numbers. 
  Usage: N50(numlist) 

  Based on the definition from this SEQanswers post 
  http://seqanswers.com/forums/showpost.php?p=7496&postcount=4 
  (modified Broad Institute's definition 
  https://www.broad.harvard.edu/crd/wiki/index.php/N50) 
   
  See SEQanswers threads for details: 
  http://seqanswers.com/forums/showthread.php?t=2857 
  http://seqanswers.com/forums/showthread.php?t=2332 
  """ 
  numlist.sort(reverse = True)
  s = sum(numlist)
  limit = s * 0.5
  for l in numlist:
      s -= l
      if s <= limit:
          return l

def count_lowercase(string):
    '''
    Tested to be quicker than "return sum(1 for c in string if c.islower())"
    
    Parameters:
        string -- any string, but typically a nucleotide sequence
    Returns:
        numLowercase -- an integer of the number of lowercase characters
                        in this string
    '''
    return sum(map(str.islower, string))

def count_char(string, char, case="upper"):
    '''
    Tested to be quicker than "return sum(1 for c in string if c.islower())"
    
    Parameters:
        string -- any string, but typically a nucleotide sequence
    Returns:
        numLowercase -- an integer of the number of lowercase characters
                        in this string
    '''
    if not case in ["upper", "lower", "any", None]:
        raise ValueError(f"num_chars doesn't understand '{case}' case value")
    
    if case == "upper":
        return string.upper().count(char)
    elif case == "lower":
        return string.lower().count(char)
    else:
        return string.count(char)
