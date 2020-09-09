import firebase from 'firebase'
import 'firebase/firestore'

// Your web app's Firebase configuration
var firebaseConfig = {
  apiKey: "AIzaSyDo0kCy6vOvX_gR9wmk5eoLBc6DMofgX2g",
  authDomain: "hsite-3e53a.firebaseapp.com",
  databaseURL: "https://hsite-3e53a.firebaseio.com",
  projectId: "hsite-3e53a",
  storageBucket: "hsite-3e53a.appspot.com",
  messagingSenderId: "1049810111052",
  appId: "1:1049810111052:web:9d93b5883ade4f954102e9"
}
// Initialize Firebase
firebase.initializeApp(firebaseConfig)

const db = firebase.firestore()
const todos = db.collection('todos')
const blog = db.collection('blog')
const auth = firebase.auth()

const increment = firebase.firestore.FieldValue.increment(1)

export default {
  db,
  auth,
  todos,
  blog,
  increment
}