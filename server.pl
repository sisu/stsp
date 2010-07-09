#!/usr/bin/perl
use IO::Socket;
use IPC::Open2;

my $sock = new IO::Socket::INET (
#               LocalHost => 'thekla',
		LocalPort => '34343',
		Proto => 'tcp',
		Listen => 1,
		Reuse => 1,
		);
print "starting loop\n";
while($cli = $sock->accept) {
	print "got conn\n";
	print $cli "lol\n";
	$s = '';
	while(<$cli>) {
		print "got data: ",$_,"\n";
		if ($_ =~ m/;/) {
			print "breaking ",$_,"\n";
			@arr = split(';', $_);
			print "arr: ",@arr," ;\n";
			$s = $s . $arr[0];
			goto out;
		}
		$s = $s.$_;
	}
out:
	my($in,$out);
	$pid = open2($in, $out, './solve', 't3.gr') || die "a.out";
	print $out ($s."\n;\n");
	while(<$in>) {
		print "res: ",$_,"\n";
		print $cli $_,"\n";
	}
	close($pid);
#       print $s."\n";
	print "closing conn\n";
	close($cli);
}
