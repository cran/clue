### <NOTE>
### It may be better to group the package-specific methods by package
### (rather than by generic).
### </NOTE>

### * Internal stuff.

.false <- function(x) FALSE
.true <- function(x) TRUE

.one_entry_per_column <-
function(x, j)
{
    ## For a matrix x and a vector of column indices j_1, ..., j_n where
    ## n is the number of rows of x, get x[1,j_1], ..., x[n,j_n].
    if(!is.matrix(x))
        stop("Argument 'x' must be a matrix.")
    x[cbind(seq(length = nrow(x)), j)]
}

".one_entry_per_column<-" <-
function(x, j, value)
{
    if(!is.matrix(x))
        stop("Argument 'x' must be a matrix.")        
    x[cbind(seq(length = nrow(x)), j)] <- value
    x
}

### * n_of_objects

## Get the number of objects in a clustering.

n_of_objects <-
function(x)
    UseMethod("n_of_objects")

### ** Default method.

n_of_objects.default <-
function(x)
    length(cl_class_ids(x))
## (Note that prior to R 2.1.0, kmeans() returned unclassed results,
## hence the best we can do for the *default* method is to look at a
## possibly existing "cluster" component.  Using the class ids incurs
## another round of method dispatch, but avoids code duplication.)

### ** Partitioning methods.

## Package stats: kmeans() (R 2.1.0 or better).
n_of_objects.kmeans <- n_of_objects.default

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
n_of_objects.partition <- n_of_objects.default

## Package cclust: cclust().
n_of_objects.cclust <- n_of_objects.default

## Package e1071: cmeans() gives objects of class "fclust".
n_of_objects.fclust <-
function(x)
    nrow(x$membership)
## Package e1071: cshell().
n_of_objects.cshell <- n_of_objects.fclust
## Package e1071: bclust().
n_of_objects.bclust <- n_of_objects.default

## Package mclust: Mclust().
n_of_objects.Mclust <- n_of_objects.default

### ** Hierarchical methods.

## Package stats: hclust().
n_of_objects.hclust <-
function(x)
    length(x$order)

## Package cluster: agnes() and diana() give objects inheriting from
## class "twins".
n_of_objects.twins <- n_of_objects.hclust
## Package cluster: mona().
n_of_objects.mona <- n_of_objects.hclust

### ** Others.

## Package clue: Ensembles.
n_of_objects.cl_ensemble <-
function(x)
    attr(x, "n_of_objects")
## Package clue: Memberships.
n_of_objects.cl_membership <- nrow
## Package clue: cl_pclust().
n_of_objects.cl_pclust <- n_of_objects.fclust
## Package clue: Ultrametrics.
n_of_objects.cl_ultrametric <-
function(x)
    attr(x, "Size")

### * n_of_classes

## Get the number of classes in a (hard or soft) partition.

## <NOTE>
## We generally allow for classes to be empty, unlike the current
## version of kmeans().  Package cclust has a version of k-means which
## does not stop in case of empty classes.
## However, we only count NON-EMPTY classes here.
## </NOTE>

n_of_classes <-
function(x)
    UseMethod("n_of_classes")

## Default method.
n_of_classes.default <-
function(x)
    length(unique(cl_class_ids(x)))

## Package stats: kmeans() (R 2.1.0 or better).
n_of_classes.kmeans <- n_of_classes.default

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
n_of_classes.fanny <-
function(x)
    sum(colSums(x$membership) > 0)
n_of_classes.partition <- n_of_classes.default

## Package cclust: cclust().
n_of_classes.cclust <- n_of_classes.default

## Package e1071: cmeans() gives objects of class "fclust".
n_of_classes.fclust <- n_of_classes.fanny
## Package e1071: cshell().
n_of_classes.cshell <- n_of_classes.fanny
## Package e1071: bclust().
n_of_classes.bclust <- n_of_classes.default

## Package mclust: Mclust().
n_of_classes.Mclust <- n_of_classes.default

## Package clue: Memberships.
n_of_classes.cl_membership <-
function(x)
    attr(x, "n_of_classes")
## Package clue: cl_pclust().
n_of_classes.cl_pclust <- n_of_classes.fanny

### * cl_class_ids

## Get ids of classes in a partition.
## <NOTE>
## Currently, all supported soft partitioning methods provide a softmax
## hard partitioning as well.
## </NOTE>

cl_class_ids <-
function(x)
    UseMethod("cl_class_ids")

## Default method.
cl_class_ids.default <-
function(x)
{
    ## Assume data structure returned by kmeans.
    x$cluster
}

## Package stats: kmeans() (R 2.1.0 or better).
cl_class_ids.kmeans <- cl_class_ids.default

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
cl_class_ids.partition <-
function(x)
    x$clustering

## Package cba: ccfkms().
cl_class_ids.ccfkms <-
function(x)
    as.vector(x$cl)
## Package cba: rockCluster() returns objects of class "rock".
cl_class_ids.rock <-
function(x)
    as.vector(x$cl)

## Package cclust: cclust().
cl_class_ids.cclust <- cl_class_ids.default

## Package e1071: cmeans() gives objects of class "fclust".
cl_class_ids.fclust <- cl_class_ids.default
## Package e1071: cshell().
cl_class_ids.cshell <- cl_class_ids.default
## Package e1071: bclust().
cl_class_ids.bclust <- cl_class_ids.default

## Package flexclust: kcca() returns objects of S4 class "kcca" which
## extends S4 class "flexclust".
cl_class_ids.kcca <-
function(x)
    flexclust::cluster(x)

## Package flexmix: class "flexmix".
cl_class_ids.flexmix <-
function(x)
    flexmix::cluster(x)

## Package kernlab: specc() and kkmeans() return objects of S4 class
## "specc".
cl_class_ids.specc <-
function(x)
    as.vector(unclass(x))

## Package mclust: Mclust().
cl_class_ids.Mclust <-
function(x)
    x$classification

## Package clue: Memberships.
cl_class_ids.cl_membership <-
function(x)
    max.col(x)
## (Cannot do cl_class_ids.cl_membership <- max.col for generic/method
## consistency.)
## Package clue: cl_pclust().
cl_class_ids.cl_pclust <- cl_class_ids.default

### * is.cl_partition

## Determine whether an object is a (generalized) partition.
## Note that this includes both hard and soft partitions, and allows
## sums of memberships of objects to be less than one.

is.cl_partition <-
function(x)
    UseMethod("is.cl_partition")

## Default method.
is.cl_partition.default <-
function(x)
{
    ## Ugly, but what else can we do?
    ## <FIXME 2.1.0>
    ## Maybe change eventually now that in R 2.1.0 or better,
    ## stats::kmeans() returns something classed ...) 
    !is.null(x$cluster)
    ## </FIXME>    
}

## Package stats: kmeans() (R 2.1.0 or better).
is.cl_partition.kmeans <- .true

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
is.cl_partition.partition <- .true

## Package cba: ccfkms().
is.cl_partition.ccfkms <- .true
## Package cba: rockCluster() returns objects of class "rock".
is.cl_partition.rock <- .true

## Package cclust: cclust().
is.cl_partition.cclust <- .true

## Package e1071: cmeans() gives objects of class "fclust".
is.cl_partition.fclust <- .true
## Package e1071: cshell().
is.cl_partition.cshell <- .true
## Package e1071: bclust().
is.cl_partition.bclust <- .true

## Package flexclust: kcca() returns objects of S4 class "kcca" which
## extends S4 class "flexclust".
is.cl_partition.kcca <- .true

## Package flexmix: class "flexmix".
is.cl_partition.flexmix <- .true

## Package kernlab: specc() and kkmeans() return objects of S4 class
## "specc".
is.cl_partition.specc <- .true

## Package mclust: Mclust().
is.cl_partition.Mclust <- .true

## Package clue: Memberships.
is.cl_partition.cl_membership <- .true
## Package clue: cl_pclust().
is.cl_partition.cl_pclust <- .true

### * as.cl_partition

## Note that cl_partition conceptually is a virtual class, so there are
## no prototypes and no cl_partition() creator.

as.cl_partition <-
function(x)
    if(is.cl_partition(x)) x else as.cl_membership(x)

### * is.cl_hard_partition

## Determine whether an object is a hard partition.

is.cl_hard_partition <-
function(x)
    UseMethod("is.cl_hard_partition")

## Default method.
is.cl_hard_partition.default <-
function(x)
{
    ## Ugly, but what else can we do?
    ## <FIXME 2.1.0>
    ## Maybe change eventually now that in R 2.1.0 or better,
    ## stats::kmeans() returns something classed ...) 
    !is.null(x$cluster)
    ## </FIXME>    
}

## Package stats: kmeans() (R 2.1.0 or better).
is.cl_hard_partition.kmeans <- .true

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
## <NOTE>
## Of course, fuzzy clustering can also give a hard partition ...
is.cl_hard_partition.fanny <-
function(x)
{
    all(rowSums(cl_membership(x)) == 1)
}
## </NOTE>
is.cl_hard_partition.partition <- .true

## Package cba: ccfkms().
is.cl_hard_partition.ccfkms <- .true
## Package cba: rockCluster() returns objects of class "rock".
is.cl_hard_partition.rock <- .true

## Package cclust: cclust().
is.cl_hard_partition.cclust <- .true

## Package e1071: cmeans() gives objects of class "fclust".
is.cl_hard_partition.fclust <- is.cl_hard_partition.fanny
## Package e1071: cshell().
is.cl_hard_partition.cshell <- is.cl_hard_partition.fanny
## Package e1071: bclust().
is.cl_hard_partition.bclust <- .true

## Package flexclust: kcca() returns objects of S4 class "kcca" which
## extends S4 class "flexclust".
is.cl_hard_partition.kcca <- .true

## Package flexmix: class "flexmix".
is.cl_hard_partition.flexmix <- is.cl_hard_partition.fanny

## Package kernlab: specc() and kkmeans() return objects of S4 class
## "specc".
is.cl_hard_partition.specc <- .true

## Package mclust: Mclust().
is.cl_hard_partition.Mclust <- is.cl_hard_partition.fanny

## Package clue: Memberships.
is.cl_hard_partition.cl_membership <-
function(x)
    attr(x, "is_cl_hard_partition")
## Package clue: cl_pclust().
is.cl_hard_partition.cl_pclust <- is.cl_hard_partition.fanny

### * as.cl_hard_partition

as.cl_hard_partition <-
function(x)
{
    if(is.cl_hard_partition(x)) x
    else if(is.cl_partition(x)) {
        ## A soft cl_partition ...
        .cl_membership_from_class_ids(labels(x)[[2]][cl_class_ids(x)])
    }
    else if(is.matrix(x)) {
        ## A matrix of raw memberships, hopefully ...
        .cl_membership_from_class_ids(labels(x)[[2]][max.col(x)])
    }
    else if(is.atomic(x)) {
        ## A vector of raw class ids, hopefully ...
        .cl_membership_from_class_ids(x)
    }
    else
        stop("Cannot coerce to 'cl_hard_partition'.")
}

### * is.cl_soft_partition

## Determine whether an object is a soft partition.

is.cl_soft_partition <-
function(x)
    is.cl_partition(x) && ! is.cl_hard_partition(x)

### * .maybe_is_proper_soft_partition

## Determine whether an object might be a proper soft partition (in the
## sense that it is a cl_partition but not a cl_hard_partition).
## This is mostly useful when computing fuzziness measures.

.maybe_is_proper_soft_partition <-
function(x)
    UseMethod(".maybe_is_proper_soft_partition")
.maybe_is_proper_soft_partition.default <- .false
.maybe_is_proper_soft_partition.fanny <- .true
.maybe_is_proper_soft_partition.fclust <- .true
.maybe_is_proper_soft_partition.cshell <- .true
.maybe_is_proper_soft_partition.flexmix <- .true
.maybe_is_proper_soft_partition.Mclust <- .true
.maybe_is_proper_soft_partition.cl_membership <-
function(x)
    !attr(x, "is_cl_hard_partition")
.maybe_is_proper_soft_partition.cl_pclust <-
function(x)
    x$m > 1

### * is.cl_hierarchy

## Determine whether an object is a hierarchy.
## Note that strictly speaking we only deal with indexed hierarchies so
## that there is always a corresponding ultrametric (and dendrogram).

is.cl_hierarchy <-
function(x)
    UseMethod("is.cl_hierarchy")

## Default method.
is.cl_hierarchy.default <- .false

## Package stats: hclust().
is.cl_hierarchy.hclust <- .true

## Package cluster: agnes() and diana() give objects inheriting from
## class "twins".
is.cl_hierarchy.twins <- .true
## Package cluster: mona().
is.cl_hierarchy.mona <- .true

## Package clue: Ultrametrics.
is.cl_hierarchy.cl_ultrametric <- .true

### * as.cl_hierarchy

## Note that cl_partition conceptually is a virtual class, so there are
## no prototypes and no cl_hierarchy() creator.

as.cl_hierarchy <-
function(x)
    if(is.cl_hierarchy(x)) x else as.cl_ultrametric(x)


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
