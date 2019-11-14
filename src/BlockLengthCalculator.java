/*
 * Script for parsing phase block information from HapCUT2 and formatting it so it can be easily plotted
 */
import java.util.*;
import java.io.*;
public class BlockLengthCalculator {
	static String HAP_FILE = "";
	static String OUT_FILE = "";
	static String CHR_LENGTHS_FILE = "";
public static void main(String[] args) throws Exception
{
	if(args.length > 0)
	{
		for(String arg : args)
		{
			int equalsIdx = arg.indexOf('=');
			if(equalsIdx == -1)
			{
				continue;
			}
			String key = arg.substring(0, equalsIdx);
			String value = arg.substring(equalsIdx+1);
			if(key.equalsIgnoreCase("hap"))
			{
				HAP_FILE = value;
			}
			else if(key.equalsIgnoreCase("out"))
			{
				OUT_FILE = value;
			}
			else if(key.equalsIgnoreCase("chrlengths"))
			{
				CHR_LENGTHS_FILE = value;
			}
		}
	}
	else
	{
		HAP_FILE = "/home/mkirsche/enc004_all.hap";
		OUT_FILE = "phase_blocks.txt";
		CHR_LENGTHS_FILE = "human_chr_lengths.txt";
	}
	if(HAP_FILE.length() == 0 || OUT_FILE.length() == 0 || CHR_LENGTHS_FILE.length() == 0)
	{
		System.out.println("Usage: java BlockLengthCalculator hap=[FILE] out=[FILE] chrlengths=[FILE]");
	}
	
	/*
	 * Get the length of each chromosome
	 */
	HashMap<String, Long> lengths = new HashMap<String, Long>();
	Scanner lengthsInput = new Scanner(new FileInputStream(new File(CHR_LENGTHS_FILE)));
	while(lengthsInput.hasNext())
	{
		String line = lengthsInput.nextLine();
		String[] tokens = line.split("\t");
		String name = "chr" + tokens[0];
		long length = Long.parseLong(tokens[2].replaceAll(",", ""));
		lengths.put(name, length);
	}
	lengthsInput.close();
	
	
	/*
	 * Read in the phase blocks and bin them by chromosome
	 */
	Scanner input = new Scanner(new FileInputStream(new File(HAP_FILE)));
	TreeMap<Name, ArrayList<Pair>> byChr = new TreeMap<Name, ArrayList<Pair>>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(!line.contains("SPAN")) continue;
		if(!input.hasNext()) break;
		
		String[] tokens = line.split(" ");
		if(tokens.length == 1) continue;
		int numVariants = Integer.parseInt(tokens[6]);
		String chr = "";
		int minPos = -1, maxPos = -1;
		for(int i = 0; i<numVariants; i++)
		{
			String after = input.nextLine();
			String[] afterTokens = after.split("\t");
			chr = afterTokens[3];
			int pos = Integer.parseInt(afterTokens[4]);
			if(minPos == -1 || pos < minPos) minPos = pos;
			if(maxPos == -1 || pos > maxPos) maxPos = pos;
		}
		
		if(!byChr.containsKey(new Name(chr)))
		{
			byChr.put(new Name(chr), new ArrayList<Pair>());
		}
		byChr.get(new Name(chr)).add(new Pair(minPos, maxPos));
	}
	
	/*
	 * For each chromosome, filter and sort the phase blocks, and print them to the file
	 */
	PrintWriter out = new PrintWriter(new File(OUT_FILE));	
	for(Name s : byChr.keySet())
	{
		if(lengths.containsKey(s.s))
		{
			System.out.println("adding end for " + s.s+" "+lengths.get(s.s));
			byChr.get(s).add(new Pair(lengths.get(s.s).intValue(), 1+lengths.get(s.s).intValue()));
		}
		if(s.s.length() > 5) continue;
		Collections.sort(byChr.get(s));
		byChr.put(s, filter(byChr.get(s)));
		out.println(s);
		for(int i = 0; i<byChr.get(s).size(); i++)
		{
			Pair val = byChr.get(s).get(i);
			out.print(val + ( i == byChr.get(s).size()-1 ? "\n" : " "));
		}
	}
	input.close();
	out.close();
}
/*
 * Sorting for chromosome names which orders by number even when the name is preceded by chr
 */
static class Name implements Comparable<Name>
{
	String s;
	Name(String ss)
	{
		s = ss;
	}
	@Override
	public int compareTo(Name o) {
		if(isNumeric() != o.isNumeric())
		{
			return isNumeric() ? -1 : 1;
		}
		if(isNumeric())
		{
			return Integer.parseInt(s.substring(3)) - Integer.parseInt(o.s.substring(3));
		}
		return s.compareTo(o.s);
	}
	
	boolean isNumeric()
	{
		for(int i = 3; i<s.length(); i++)
		{
			char c = s.charAt(i);
			if(c > '9' || c < '0') return false;
		}
		return true;
	}
	
	public String toString()
	{
		return s;
	}
}
/*
 * Remove blocks which are entirely contained in other ones
 */
static ArrayList<Pair> filter(ArrayList<Pair> pairs)
{
	ArrayList<Pair> res = new ArrayList<Pair>();
	int maxEnd = -1;
	for(int i = 0; i<pairs.size(); i++)
	{
		if(pairs.get(i).b <= maxEnd) continue;
		res.add(pairs.get(i));
		maxEnd = pairs.get(i).b;
	}
	return res;
}
static class Pair implements Comparable<Pair>
{
	int a, b;
	Pair(int aa, int bb)
	{
		a = aa; b = bb;
	}
	@Override
	public int compareTo(Pair o) {
		if(a != o.a) return a - o.a;
		return b - o.b;
	}
	public String toString()
	{
		return a + " " + b;
	}
}
}
