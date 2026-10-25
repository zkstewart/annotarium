# Copyright (C) 2026 Zachary Kenneth Stewart

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
