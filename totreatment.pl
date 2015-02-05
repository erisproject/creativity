use strict;
use warnings;

chomp(my $header = <>);
my @fields = split /,/, $header;
if (@fields != 69) { die "Wrong file format"; }

my @y;
my ($first_pre, $last_pre, $first_post, $last_post, $last_bp);
for my $i (0 .. $#fields) {
    if ($fields[$i] =~ /^pre_/) {
        $last_pre = $i;
        $first_pre = $i if not defined $first_pre;
        push @y, $fields[$i] =~ s/^pre_//r;
    }
    elsif ($fields[$i] =~ /^post_/) {
        if ($fields[$i] eq 'post_books_pirated') { $last_bp = $i; }
        $last_post = $i;
        $first_post = $i if not defined $first_post;
    }
}

if (not defined $first_pre or not defined $first_post) { die "Couldn't find first pre_/post_"; }
if ($last_pre - $first_pre != $last_post - $first_post - 1) { die "#pre != #post-1"; }
if ($last_pre >= $first_post) { die "pre after post"; }
if (not defined $last_bp) { die "didn't find post_books_pirated"; }

my @keep_fields = @fields[0 .. $first_pre-1];
for (1 .. $first_pre-1) { # (Don't dupe 0: it's the source filename)
    push @keep_fields, "pX$fields[$_]";
}
push @keep_fields, "piracy";
push @keep_fields, @y;

print join ",", @keep_fields;
print "\n";

while (<>) {
    my @f = split /,/;
    my @d1 = (@f[0 .. $first_pre-1], (0) x ($first_pre), @f[$first_pre .. $last_pre]);
    my @d2 = ($f[0], (@f[1 .. $first_pre-1]) x 2, 1, @f[$first_post .. $last_bp-1, $last_bp+1 .. $last_post]);
    print join ",", @d1;
    print "\n";
    print join ",", @d2;
    print "\n";
}
