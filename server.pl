#!/usr/bin/perl
use IO::Socket;
use IPC::Open2;
use threads;

sub handle {
	my ($cli) = @_;
	$cli->autoflush(1);
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
	$pid = open2($in, $out, './solve', 't3-2.gr') || die "a.out";
	print $out ($s."\n;\n");
	while(<$in>) {
		print "res: ",$_;
		print $cli $_;
	}
	close($pid);
#       print $s."\n";
	print "closing conn\n";
	close($cli);
}

my $sock = new IO::Socket::INET (
#               LocalHost => 'thekla',
		LocalPort => '34343',
		Proto => 'tcp',
		Listen => 1,
		Reuse => 1,
		);
print "starting loop\n";
while($cli = $sock->accept) {
#	handle($cli);
#	$t = threads->new(\&handle, $cli);
	if (($c = fork()) != 0) {
		handle($cli);
		kill("TERM", $c);
	}
	print "asdasd\n";
}
close($sock);
