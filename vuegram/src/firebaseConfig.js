import firebase from 'firebase'
import 'firebase/firestore'

const src="https://www.gstatic.com/firebasejs/7.15.4/firebase-app.js"
const config = {
    apiKey: "AIzaSyAqPfxKOU0JIT9MXRkeH0Goqh6YtexyI8E",
    authDomain: "vuegram-81295.firebaseapp.com",
    databaseURL: "https://vuegram-81295.firebaseio.com",
    projectId: "vuegram-81295",
    storageBucket: "vuegram-81295.appspot.com",
    messagingSenderId: "1092715729952",
    appId: "1:1092715729952:web:713835c4cd77cf7f349845"
}
firebase.initializeApp(config)

const db = firebase.firestore()
const auth = firebase.auth()
const currentUser = auth.currentUser

const settings = {
}
db.settings(settings)

const usersCollection = db.collection('users')
const postsCollection = db.collection('posts')
const commentsCollection = db.collection('comments')
const likesCollection = db.collection('likes')

export {
	db,
	auth,
	currentUser,
	usersCollection,
	postsCollection,
	commentsCollection,
	likesCollection
}