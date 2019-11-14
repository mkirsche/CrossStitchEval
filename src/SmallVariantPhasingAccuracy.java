/*
 * Small script for measuring the accuracy of small variant calling/phasing
 */
import java.util.*;
import java.io.*;
public class SmallVariantPhasingAccuracy {
	static String TRUTH_FILE="";
	static String PHASED_VARIANTS="";
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
			if(key.equalsIgnoreCase("truth"))
			{
				TRUTH_FILE = value;
			}
			else if(key.equalsIgnoreCase("phased_variants"))
			{
				PHASED_VARIANTS = value;
			}
		}
	}
	else
	{
		TRUTH_FILE = "gold.vcf";
		PHASED_VARIANTS = "mine.vcf";
	}
	if(TRUTH_FILE.length() == 0 || PHASED_VARIANTS.length() == 0)
	{
		System.out.println("Usage: java SmallVariantPhasingAccuracy truth=[FILE] phased_variants=[FILE]");
	}
	
	/*
	 * Read variants from ground truth file
	 */
	Scanner input = new Scanner(new FileInputStream(new File(TRUTH_FILE)));
	ArrayList<Variant> trueVariants = new ArrayList<Variant>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0 || line.startsWith("#"))
		{
			continue;
		}
		Variant var = new Variant(line);
		trueVariants.add(var);
	}
	
	System.err.println("Found " + trueVariants.size() + " variants in the ground truth");
	
	/*
	 * Read variants from callset
	 */
	ArrayList<Variant> calledVariants = new ArrayList<Variant>();
	input = new Scanner(new FileInputStream(new File(PHASED_VARIANTS)));
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0 || line.startsWith("#"))
		{
			continue;
		}
		Variant var = new Variant(line);
		calledVariants.add(var);
	}
	
	System.err.println("Found " + calledVariants.size() + " variants in the callset");
	
	/*
	 * Sort both datasets so they can be compared quickly
	 */
	Collections.sort(trueVariants);
	Collections.sort(calledVariants);
	
	
	/*
	 * Scan through both files at once, keeping track of the current point in each file
	 */
	int truthIdx = 0;
	int calledIdx = 0;
	
	TreeMap<String, ChromosomeData> cd = new TreeMap<String, ChromosomeData>();
	
	while(truthIdx < trueVariants.size() || calledIdx < calledVariants.size())
	{
		// Case 1 - already reached the end of the ground truth so false positive
		if(truthIdx == trueVariants.size())
		{
			String name = calledVariants.get(calledIdx).chr;
			if(!cd.containsKey(name))
			{
				cd.put(name, new ChromosomeData());
			}
			cd.get(name).falsePositives++;
			calledIdx++;
		}
		// Case 2 - already reached the end of the callset, so false negative
		else if(calledIdx == calledVariants.size())
		{
			String name = trueVariants.get(truthIdx).chr;
			if(!cd.containsKey(name))
			{
				cd.put(name, new ChromosomeData());
			}
			cd.get(name).missed++;
			cd.get(name).missedByGt[trueVariants.get(truthIdx).genotype + 1]++;
			cd.get(name).count++;
			truthIdx++;
		}
		else
		{
			Variant called = calledVariants.get(calledIdx);
			Variant truth = trueVariants.get(truthIdx);
			// Case 3 - called variant is at earlier position so has no match and is false positive
			if(called.compareTo(truth) < 0)
			{
				String name = calledVariants.get(calledIdx).chr;
				if(!cd.containsKey(name))
				{
					cd.put(name, new ChromosomeData());
				}
				cd.get(name).falsePositives++;
				calledIdx++;
			}
			// Case 4 - true variant is at earlier position so has no match and is false negative
			else if(called.compareTo(truth) > 0)
			{
				String name = trueVariants.get(truthIdx).chr;
				if(!cd.containsKey(name))
				{
					cd.put(name, new ChromosomeData());
				}
				cd.get(name).missed++;
				cd.get(name).missedByGt[trueVariants.get(truthIdx).genotype + 1]++;
				cd.get(name).count++;
				truthIdx++;
			}
			// Case 5 - match, so compare ref/alt/genotype info
			else
			{
				String name = trueVariants.get(truthIdx).chr;
				if(!cd.containsKey(name))
				{
					cd.put(name, new ChromosomeData());
				}
				if(!called.ref.equalsIgnoreCase(truth.ref) || !called.alt.equalsIgnoreCase(truth.alt))
				{
					cd.get(name).wrongRefAlt++;
					cd.get(name).missedByGt[trueVariants.get(truthIdx).genotype + 1]++;
				}
				else
				{
					if(truth.genotype != called.genotype)
					{
						cd.get(name).wrongPhasing++;
					}
					cd.get(name).gtCalls[truth.genotype+1][called.genotype+1]++;
					if(truth.genotype == 1 && called.genotype == 2)
					{
						cd.get(name).phaseString.append("0");
					}
					else if(truth.genotype == 2 && called.genotype == 1)
					{
						cd.get(name).phaseString.append("0");
					}
					else if(truth.genotype == 1 && called.genotype == 1)
					{
						cd.get(name).phaseString.append("1");
					}
					else if(truth.genotype == 2 && called.genotype == 2)
					{
						cd.get(name).phaseString.append("1");
					}
				}
				cd.get(name).count++;
				truthIdx++;
				calledIdx++;
			}
		}
	}
	
	for(String s : cd.keySet())
	{
		System.err.println("Results for " + s);
		ChromosomeData cur = cd.get(s);
		System.err.println("Total of " + cur.count + " variants in the ground truth");
		System.err.println("Failed to call " + cur.missed + " variants in the ground truth");
		System.err.println("Wrong ref and/or alt for " + cur.wrongRefAlt + " variants");
		System.err.println("Wrong genotype for " + cur.wrongPhasing + " variants");
		System.err.println("Called " + cur.falsePositives + " false positives");
		System.err.println("Genotype confusion matrix (ungenotyped, 0|0, 0|1, 1|0, 1|1 - row is truth and column is call)");
		for(int i = 0; i<cur.gtCalls.length; i++)
		{
			for(int j = 0; j<cur.gtCalls[i].length; j++)
			{
				System.err.print(cur.gtCalls[i][j] + ((j == cur.gtCalls[i].length - 1) ? "\n" : " "));
			}
		}
		int[] ser = cur.switchErrorRate();
		System.err.println("Switch error rate: " + 1.0 * ser[0] / ser[1]);
		System.err.println();
	}
	PrintWriter err = new PrintWriter(System.err);
	getTotalMetrics(cd, err);
	input.close();
	err.close();
}
/*
 * Get genome-wide metrics from per-chromosome metrics
 */
static void getTotalMetrics(TreeMap<String, ChromosomeData> cd, PrintWriter out)
{
	int[] totalSer = new int[2];
	int[][] totalGtCalls = new int[5][5];
	int totalCount = 0;
	int totalFp = 0;
	int totalFn = 0;
	int totalWrongRefAlt = 0;
	int[] totalMissedByGt = new int[totalGtCalls.length];
	for(String s : cd.keySet())
	{
		ChromosomeData cur = cd.get(s);
		int[] ser = cur.switchErrorRate();
		totalSer[0] += ser[0];
		totalSer[1] += ser[1];
		for(int i = 0; i<totalGtCalls.length; i++)
		{
			for(int j = 0; j<totalGtCalls[i].length; j++)
			{
				totalGtCalls[i][j] += cur.gtCalls[i][j];
			}
		}
		for(int i = 0; i<totalMissedByGt.length; i++)
		{
			totalMissedByGt[i] += cur.missedByGt[i];
		}
		totalCount += cur.count;
		totalFp += cur.falsePositives;
		totalFn += cur.missed;
		totalWrongRefAlt += cur.wrongRefAlt;
	}
	int[] rowTots = new int[totalGtCalls.length];
	int[] colTots = new int[totalGtCalls[0].length];
	for(int i = 0; i<totalGtCalls.length; i++)
		for(int j = 0; j<totalGtCalls[i].length; j++)
		{
			rowTots[i] += totalGtCalls[i][j];
			colTots[j] += totalGtCalls[i][j];
		}
	out.println();
	out.println("Detection accuracy:");
	out.println("  Total count in ground truth: " + totalCount);
	out.println("  Number of false positives: " + totalFp);
	out.println("  Number of false negatives: " + totalFn);
	out.println("  Number of matches with wrong ref/alt: " + totalWrongRefAlt);
	out.println();
	out.println("Phasing Accuracy");
	out.println("  Homozygous variants in ground truth: " + (rowTots[4] + totalMissedByGt[4]));
	out.println("    Called homozygous: " + totalGtCalls[4][4]);
	out.println("    Called as phased heterozygous: " + (totalGtCalls[4][2] + totalGtCalls[4][3]));
	out.println("    Called as unphased heterozygous: " + totalGtCalls[4][0]);
	out.println("    Called as 0|0: " + totalGtCalls[4][1]);
	out.println("    Not called: " + totalMissedByGt[4]);
	out.println();
	out.println("  Heterozygous variants in ground truth: " + (rowTots[2] + rowTots[3] + totalMissedByGt[2] + totalMissedByGt[3]));
	out.println("    Called homozygous: " + (totalGtCalls[2][4] + totalGtCalls[3][4]));
	out.println("    Called as phased heterozygous: " 
	    + (totalGtCalls[2][2] + totalGtCalls[2][3] + totalGtCalls[3][2] + totalGtCalls[3][3]));
	out.println("    Called as unphased heterozygous: " + (totalGtCalls[2][0] + totalGtCalls[3][0]));
	out.println("    Called as 0|0: " + (totalGtCalls[2][1] + totalGtCalls[3][1]));
	out.println("    Not called: " + (totalMissedByGt[2] + totalMissedByGt[3]));
	out.println();
	out.println("  Switch Error Rate: " + totalSer[0] + " out of " + totalSer[1] + " = " + 1.0 * totalSer[0] / totalSer[1]);
}
/*
 * Metrics for phasing performance along a single chromosome
 */
static class ChromosomeData
{
	int falsePositives;
	int wrongRefAlt;
	int wrongPhasing;
	int[][] gtCalls; // Confusion matrix for genotype calls - rows are truth and cols are calls
	int missed;
	int count;
	StringBuilder phaseString;
	int[] missedByGt; // How many true calls of each genotype were missed
	ChromosomeData()
	{
		falsePositives = 0;
		missed = 0;
		count = 0;
		wrongRefAlt = 0;
		wrongPhasing = 0;
		phaseString = new StringBuilder("");
		gtCalls = new int[5][5];
		missedByGt = new int[gtCalls.length];
	}
	int[] switchErrorRate()
	{
		int total = 0;
		int bad = 0;
		String s = phaseString.toString();
		for(int i = 0; i<s.length()-1; i++)
		{
			char cur = s.charAt(i);
			char next = s.charAt(i+1);
			total++;
			if(cur != next) bad++;
		}
		return new int[] {bad, total};
	}
}
/*
 * Info for a small variant call
 */
static class Variant implements Comparable<Variant>
{
	String chr;
	int pos;
	String ref;
	String alt;
	int genotype;
	Variant(String line)
	{
		String[] tokens = line.split("\t");
		chr = tokens[0];
		pos = Integer.parseInt(tokens[1]);
		ref = tokens[3];
		alt = tokens[4];
		genotype = -1;
		if(tokens.length >= 10)
		{
			String gtString = tokens[9];
			int colonIdx = gtString.indexOf(':');
			if(colonIdx != -1)
			{
				genotype = gtStringToInt(gtString.substring(0, colonIdx));
			}
		}
	}
	static int gtStringToInt(String s)
	{
		if(s.equals("1|1") || s.equals("1/1"))
		{
			return 3;
		}
		if(s.equals("1|0"))
		{
			return 2;
		}
		if(s.equals("0|1"))
		{
			return 1;
		}
		if(s.equals("1/0") || s.equals("0/1"))
		{
			return -1;
		}
		return 0;
	}
	@Override
	public int compareTo(Variant o) {
		if(!chr.equals(o.chr)) return chr.compareTo(o.chr);
		return pos - o.pos;
	}
}
}
