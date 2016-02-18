;;;; distributions.asd

(asdf:defsystem #:distributions
  :description "Describe distributions here"
  :author "Your Name <your.name@example.com>"
  :license "Specify license here"
  :depends-on (:cl-quad :cl-numerics-utils)
  :serial t
  :components ((:file "package")
               (:file "distributions")))

