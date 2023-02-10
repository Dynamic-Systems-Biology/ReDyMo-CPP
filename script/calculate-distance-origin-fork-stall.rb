# Script to generate the part of the comparidon table of distances between origin and fork stall
gem 'parallel'
gem 'descriptive_statistics'
require 'parallel'
require 'descriptive_statistics'

results_dir = File.expand_path(ARGV[0])
n_rounds = 10
n_forks = (10..100).step(5)
periods = [
  0,
  100,
  1000,
  10_000,
  100_000
]

def streak_size(line)
  if line.include? '-' and line.include? 'x'
    streak = line.split('-').map(&:to_i)
    repeat = line.split('x')[1]
    kind = 'both'
  elsif line.include? '-'
    streak = line.split('-').map(&:to_i)
    kind = 'sequence'
  elsif line.include? 'x'
    repeat = line.split('x')[1]
    kind = 'repeat'
  end

  streak_len = 0
  case kind
  when 'sequence'
    streak_len = (streak[0] - streak[1]).abs
    streak_len
  when 'repeat'
  when 'both'
  end
end

p 'number of forks, steps per iteration, total iterations, average inter origin distance'
Parallel.each(n_forks) do |f|
  Parallel.each(periods) do |period|
    p "Forks: #{f} | Period: #{period}"
    streaks = []
    (0..n_rounds - 1).each do |round|
      (0..1000 - 1).each do |sim|
        ('01'..'11').each do |chromosome|
          chrm_file = "#{results_dir}/round_#{round}_false_#{f}_#{period}/simulation_#{sim}/Tb927_#{chromosome}_v5.1.cseq"
          File.open(chrm_file).each_line do |line|
            len = streak_size(line)
            streaks.push(len) if len
          end
        rescue Errno::ENOENT => e
        # p "#{chrm_file} not found"
        rescue e
          p e
        end
      end
    end
    p "#{f}, #{period}, #{streaks.median}, #{streaks.mean}, #{streaks.standard_deviation}"
  end
end
