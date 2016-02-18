;;;; package.lisp

(defpackage #:distributions
  (:use #:cl #:cl-quad #:cl-numerics-utils)
  (:export
   #:cumulative-distribution #:distribution-density
   #:distribution-mean #:distribution-variance
   #:distribution #:lower-bound #:upper-bound #:mean #:variance
   #:ggs #:top-size #:alpha
   #:rosin-rammler #:size-63.2%
   #:truncation #:x-max
   #:truncation-forward #:truncation-reverse #:truncation-reverse-deriv
   #:truncated-distribution
   #:truncated-rosin-rammler))

