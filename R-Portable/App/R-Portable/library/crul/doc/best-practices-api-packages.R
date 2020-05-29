## ----eval=FALSE----------------------------------------------------------
#  xGET <- function(url, path, args = list(), ...) {
#    cli <- crul::HttpClient$new(url, opts = list(...))
#    res <- cli$get(path = path, query = args)
#    res$raise_for_status()
#    res$raise_for_ct_json()
#    res$parse("UTF-8")
#  }

## ----eval=FALSE----------------------------------------------------------
#  x <- xGET("https://httpbin.org", "get", args = list(foo = "bar"))
#  # parse the JSON to a list
#  jsonlite::fromJSON(x)
#  # more parsing

## ----eval=FALSE----------------------------------------------------------
#  xGET("https://xxx.org", args = list(foo = "bar"), verbose = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  con <- HttpClient$new("https://httpbin.org/status/404")
#  res <- con$get()

## ----eval=FALSE----------------------------------------------------------
#  x <- fauxpas::find_error_class(res$status_code)$new()
#  #> <HTTPNotFound>
#  #>  behavior: stop
#  #>  message_template: {{reason}} (HTTP {{status}})
#  #>  message_template_verbose: {{reason}} (HTTP {{status}}).\n - {{message}}

## ----eval=FALSE----------------------------------------------------------
#  x$do(res)
#  #> Error: Not Found (HTTP 404)

## ----eval=FALSE----------------------------------------------------------
#  x$do(res, template = "{{status}}\n  --> {{reason}}")
#  #> Error: 404
#  #>  --> Not Found

## ----eval=FALSE----------------------------------------------------------
#  x$do_verbose(res)
#  #> Error: Not Found (HTTP 404).
#  #>  - The server has not found anything matching the Request-URI. No indication
#  #>  is given of whether the condition is temporary or permanent. The 410 (Gone)
#  #>  status code SHOULD be used if the server knows, through some internally configurable
#  #>  mechanism, that an old resource is permanently unavailable and has no forwarding
#  #>  address. This status code is commonly used when the server does not wish to
#  #>  reveal exactly why the request has been refused, or when no other response
#  #> is applicable.

## ----eval=FALSE----------------------------------------------------------
#  x$behavior <- "warning"
#  x$do(res)
#  #> Warning message:
#  #> Not Found (HTTP 404)
#  x$behavior <- "message"
#  x$do(res)
#  #> Not Found (HTTP 404)

