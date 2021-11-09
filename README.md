Here's an FFT implementation I used many years ago to implement a spectrum analyzer
in a music recorder. It worked great back then, should work even better now!

Found the original at: https://www.musicdsp.org/en/latest/Analysis/ in FFT section:

--->

A paper (postscript) and some C++ source for 4 different fft algorithms, compiled by Toth
Laszlo from the Hungarian Academy of Sciences Research Group on Artificial Intelligence.
Toth says: "I've found that Sorensen's split-radix algorithm was the fastest, so I use
this since then (this means that you may as well delete the other routines in my source -
if you believe my results)."

<---

So I put those 150 lines directly into my recorder and forgot the rest.
This is the code I mentioned at NDC TechTown 2021 in my talk about 'Real Programming'.
Included in a separate file for easy access. 

I just love the original coding style, the haphazard formatting, the lack of comments,
it's everything that good code should be - fast, consice and without formal crap. Also
not really seaplusplus, this is plain old C the way we like it.


Sjur Julin

