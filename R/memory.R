memory <- utils::hashtab()

prepareMemory <- function() {
  utils::clrhash(memory)
}

addToMemory <- function(key, value) {
  utils::sethash(memory, key, value)
}

getFromMemory <- function(key) {
  utils::gethash(memory, key, nomatch = NULL)
}
