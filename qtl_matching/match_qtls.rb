$interestingQTL = ARGV[0]
$suspectGenes = ARGV[1]
$output = ARGV[2]

$regExForPosition = /([A-Z0-9]*):([0-9]*)-([0-9]*)/

$exceptions = {
  "LOC100694761" => "LWS",
  "LOC100695018" => "SWS2B",
  "LOC100695287" => "SWS2A",
  "LOC100710676" => "RH2A-beta",
  "LOC100710942" => "RH2A-alpha",
  "LOC100711209" => "RH2B",
  "LOC100710249" => "SWS1",
  "LOC100711443" => "rho",
  "LOC100702555" => "Tbx2a",
  "LOC100705124" => "Rx1",
};

def main()
  interesting_genes = {}
  gene_info = {}
  outPointer = File.open($output,"w")
  outPointer.puts("Modules of interest");
  puts($suspectGenes)
  File.open($suspectGenes).each_line{|line|
    arr = line.split("\t")
    if(arr.length > 4 && arr[1] =~ $regExForPosition)
      interesting_genes[arr[0]] =Locus.new($1,$2,$3)
      gene_info[arr[0]] = arr.join("\t")
      if($exceptions.keys.include?(arr[0]))
        outPointer.puts(line.chomp() + "\t" + $exceptions[arr[0]]);
      end
    end
  }
  outPointer.puts();
  File.open($interestingQTL,"r").each_line{|line|
    line.chomp!
    loci = line.split("\t")
    gene = loci.shift
    loci.each{|locus|
      if(locus =~ $regExForPosition)
        comparison = Locus.new($1,$2,$3)
        outPointer.puts(">" +gene + ":" + comparison.toString())
        interesting_genes.each{|geneOut,geneLocus|
          if(geneLocus.overlaps?(comparison))
            outPointer.puts(gene_info[geneOut])
          end
        }
        outPointer.puts("");
      end
    }
  }
  outPointer.close
end
class Locus

  attr_accessor :linkageGroup
  attr_accessor :startPos
  attr_accessor :endPos
  
  def initialize(lg, start, fin)
    @linkageGroup = lg
    @startPos = start.to_i()
    @endPos = fin.to_i()
  end

  def overlaps?(otherLocus)
    return (otherLocus.linkageGroup == @linkageGroup && !(otherLocus.startPos > @endPos || otherLocus.endPos < @startPos))
  end

  def toString()
    return "#{@linkageGroup}:#{@startPos}-#{@endPos}";
  end
  
end

main();
