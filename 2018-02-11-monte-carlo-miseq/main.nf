
params.simulator = '/home/chm2059/chm2059/1216_bishoy/results/experiments/040617_monteCarlo/main.py'

simulations = Channel
  .from(
    [278,13,0.00025],
    [278,14,0.0005],
    [278,15,0.001],
    [1392,11,0.00025],
    [1392,12,0.0005],
    [1392,13,0.001],
    [2785,10,0.00025],
    [2785,11,0.0005],
    [2785,12,0.001],
    [13927,8,0.00025],
    [13927,9,0.0005],
    [13927,10,0.001],
    [27855,7,0.00025],
    [27855,8,0.0005],
    [27855,9,0.001],
    [65710,5,0.00025],
    [65710,6,0.0005],
    [65710,7,0.001]
  )
  .spread([0.8,0.9])
  .spread(["1e-5","1e-6","1e-7"])

process simulate {

  executor 'sge'
  clusterOptions '-l h_vmem=8G -pe smp 10 -l os=rhel6.3'
  storeDir "${baseDir}/simulation/"

  input:
    set inputAmt, repliCycles, sampling, pcrEffic, polymeraseError from simulations

  output:
    file("${inputAmt}_${polymeraseError}_${polymeraseError}_${repliCycles}_20_${pcrEffic}_${pcrEffic}_${sampling}.csv") into simulate_out

  """
  python ${params.simulator} \
    --initpop ${inputAmt} \
    --polymerase_error_truseq ${polymeraseError} \
    --polymerase_error_replig ${polymeraseError} \
    --cycles_replig ${repliCycles} \
    --cycles_truseq 18 \
    --pcr_efficiency_truseq ${pcrEffic} \
    --pcr_efficiency_replig ${pcrEffic} \
    --replig_sampling ${sampling} \
    --prefix ${inputAmt}_${polymeraseError}_${polymeraseError}_${repliCycles}_20_${pcrEffic}_${pcrEffic}_${sampling} \
    --p 10 \
    --n 100

  """
}
