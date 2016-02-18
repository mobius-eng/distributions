(in-package #:distributions)

;; ** Distributions protocol
(defgeneric cumulative-distribution (distribution point))
(defgeneric distribution-density (distribution point))
(defgeneric distribution-mean (distribution))
(defgeneric distribution-variance (distribution))

;; ** Generic distribution
(defclass distribution ()
  ((lower-bound
    :initform +double-float-minus-infinity+
    :reader distribution-lower-bound
    :documentation
    "Lower bound of a distribution,
     defaults to +DOUBLE-FLOAT-MINUS-INFINITY+")
   (upper-bound
    :initform +double-float-plus-infinity+
    :reader distribution-upper-bound
    :documentation
    "Upper bound of a distribution,
     defaults to +DOUBLE-FLOAT-PLUS-INFINITY+")
   (mean
    :documentation
    "Distribution mean. If left unbound, the first call to
     DISTRIBUTION-MEAN will calculate it numerically and bind
     the slot to the result")
   (variance
    :documentation
    "Distribution vairance. If left unbound, the first call to
     DISTRIBUTION-MEAN will calculate it numerically and bind
     the slot to the result"))
  (:documentation
   "Basic distribution structure providing the following slots:
    LOWER-BOUND : lower bound of a distribution,
                  defaults to +DOUBLE-FLOAT-MINUS-INFINITY+
    UPPER-BOUND : upper bound of a distribution,
                  defaults to +DOUBLE-FLOAT-PLUS-INFINITY+
    MEAN        : distribution mean. If left unbound, the first call to
                  DISTRIBUTION-MEAN will calculate it numerically and bind
                  the slot to the result
    VARIANCE    : distribution vairance. If left unbound, the first call to
                  DISTRIBUTION-MEAN will calculate it numerically and bind
                  the slot to the result"))

(defmethod distribution-mean ((dist distribution))
  "Default numerical calculation of the mean. If MEAN slot is bound,
   just returns it. Requires the distribution to have implemented
   DISTRIBUTION-DENSITY"
  (if (slot-boundp dist 'mean)
      (slot-value dist 'mean)
      (with-slots (lower-bound upper-bound mean) dist
        (setf mean (integrate (lambda (x)
                                (* x (distribution-density dist x)))
                              lower-bound
                              upper-bound)))))


(defmethod distribution-variance ((dist distribution))
  "Default numerical calculation of the variabce. If VARIANCE slot is bound,
   just returns it. Requires the distribution to have implemented
   DISTRIBUTION-DENSITY"
  (if (slot-boundp dist 'variance)
      (slot-value dist 'variance)
      (with-slots (lower-bound upper-bound variance) dist
        (let ((mean (distribution-mean dist)))
          (setf variance (integrate (lambda (x)
                                      (* (expt (- x mean) 2)
                                         (distribution-density dist x)))
                                    lower-bound
                                    upper-bound))))))

;; ** GGS
(defclass ggs (distribution)
  ((top-size :initarg :top-size :reader top-size)
   (alpha :initarg :alpha :reader alpha))
  (:documentation
   "Gates-Gaudin-Schuhmann distribution:

                 a
           / x \\ 
    P(x) = | - |
           \\ x'/

"))

(defmethod initialize-instance :after ((obj ggs) &key)
  (with-slots (lower-bound upper-bound mean variance top-size alpha) obj
    (setf lower-bound 0d0)
    (setf upper-bound top-size)
    (setf mean (/ (* alpha top-size) (1+ alpha)))
    (setf variance (/ (* alpha (expt top-size 2))
                      (* (expt (1+ alpha) 2) (+ alpha 2))))))

(defun ggs (top-size alpha)
  "Constructs Gates-Gaudin-Schuhmann with given parameters"
  (make-instance 'ggs :top-size top-size :alpha alpha))

(defmethod cumulative-distribution ((distribution ggs) point)
  (cond ((minusp point) 0)
        ((> point (top-size distribution)) 1)
        (t (expt (/ point (top-size distribution)) (alpha distribution)))))

(defmethod distribution-density ((distribution ggs) point)
  (cond ((or (minusp point) (> point (top-size distribution))) 0)
        (t (with-slots (top-size alpha) distribution
             (* alpha
                (expt top-size (- alpha))
                (expt (/ point top-size) (1- alpha)))))))

;; ** Rosin-Rammler
(defclass rosin-rammler (distribution)
  ((size-63.2% :initarg :size-63.2% :reader size-63.2%)
   (alpha :initarg :alpha :reader alpha))
  (:documentation
   "Rosin-Rammler (Weibull) distribution:

                   _                _
                  |  /          \\ a  |
                  |  |     x    |    |
   P(x) = 1 - exp |- | -------- |    |
                  |  |   x      |    |
                  |_ \\    63.2  /   _|
                   
"))

(defmethod initialize-instance :after ((obj rosin-rammler) &key)
  (with-slots (lower-bound alpha size-63.2% mean variance) obj
    (setf lower-bound 0d0)))

(defmethod cumulative-distribution ((dist rosin-rammler) point)
  (cond ((minusp point) 0)
        (t (with-slots (size-63.2% alpha) dist
             (- 1 (exp (- (expt (/ point size-63.2%) alpha))))))))

(defmethod distribution-density ((dist rosin-rammler) point)
  (with-slots (alpha size-63.2%) dist
    (let* ((fraction (/ point size-63.2%))
           (arg (expt fraction (1- alpha))))
      (* (exp (- (* fraction arg))) alpha arg (/ size-63.2%)))))

(defun rosin-rammler (size-63.2% alpha)
  "Construct Rosin-Rammler (Weibull) distribution with given parameters"
  (make-instance 'rosin-rammler :size-63.2% size-63.2% :alpha alpha))

;; ** Truncation

(defclass truncation ()
  ((x-max :initarg :x-max
          :reader x-max
          :documentation "Max value of the truncated quantity"))
  (:documentation
   "Top value truncation"))

(defun truncation-forward (truncation unbound-x)
  "Produces bound (0, x-max) value from unbound"
  (* (x-max truncation) (/ unbound-x (+ (x-max truncation) unbound-x))))

(defun truncation-reverse (truncation bound-y)
  "Produces unbound value from bound"
  (with-slots (x-max) truncation
    (* x-max (/ bound-y (- x-max bound-y)))))

(defun truncation-reverse-deriv (truncation bound-y)
  "Derivative of the reverse truncation function"
  (with-slots (x-max) truncation
    (expt (/ x-max (- x-max bound-y)) 2)))

(defun truncation (x-max)
  "Constructs truncation"
  (make-instance 'truncation :x-max x-max))

;; ** Truncated Rosin-Rammler
(defclass truncated-distribution (distribution)
  ((distribution :initarg :distribution)
   (truncation :initarg :truncation))
  (:documentation
   "Truncated distribution with capped top size"))

(defmethod initialize-instance :after ((obj truncated-distribution) &key)
  (with-slots (lower-bound upper-bound truncation) obj
    (setf lower-bound 0d0)
    (setf upper-bound (x-max truncation))))

(defmethod distribution-density ((dist truncated-distribution) point)
  (with-slots (distribution truncation) dist
    (* (distribution-density distribution (truncation-reverse truncation point))
       (truncation-reverse-deriv truncation point))))


(defmethod cumulative-distribution ((dist truncated-distribution) point)
  (with-slots (distribution truncation) dist
    (cumulative-distribution distribution (truncation-reverse truncation point))))

(defun truncated-rosin-rammler (size-63.2% alpha x-max)
  "Constructs truncated Rosin-Rammler distribution"
  (let ((x-63.2 (/ (* size-63.2% x-max) (- x-max size-63.2%))))
   (make-instance 'truncated-distribution
       :distribution (rosin-rammler x-63.2 alpha)
       :truncation (truncation x-max))))
